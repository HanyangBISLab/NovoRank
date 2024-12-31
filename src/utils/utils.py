import os
import logging
import numpy as np
import pandas as pd

from tqdm import tqdm
from sklearn.cluster import DBSCAN
from collections import defaultdict

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def top1(dataset) -> pd.DataFrame:

    """
    This ensures only the top-ranked result among the top N candidates is retained.

    Args:
        dataset (pd.DataFrame): Input dataset containing ranked results.

    Returns:
        pd.DataFrame: Dataset containing only the top 1 ranked result.
    """

    return dataset.drop_duplicates(subset=['Source File', 'Scan number'], keep='first').reset_index(drop=True)


def drop_missing_cluster_info(dataset) -> pd.DataFrame:

    """
    Removes rows from the dataset where clustering information is missing.

    Args:
        dataset (pd.DataFrame): Input dataset.

    Returns:
        pd.DataFrame: Dataset with rows containing clustering information.
    """

    return dataset.dropna(subset=['Cluster']).reset_index(drop=True)


def labeling(dataset) -> pd.DataFrame:

    """
    Adds a 'Label' column to the dataset by comparing the 'Sequence' and 'GT' (Ground Truth) columns.

    Args:
        dataset (pd.DataFrame): A DataFrame containing 'Sequence' and 'GT' columns.

    Returns:
        pd.DataFrame: The original DataFrame with an added 'Label' column.

    Labeling logic:
        - 0 if 'Sequence' and 'GT' are the same (match).
        - 1 if 'Sequence' and 'GT' are different (mismatch).
    """

    dataset['Label'] = (dataset['Sequence'] != dataset['GT']).astype(int)

    return dataset


class ClusterRefinement:
    def __init__(self, cluster_info, mgf_info):

        """
        Initialize the ClusterRefinement class with a dataset.
        Merges information and cluster data into a single dataset.
        """

        self.dataset = pd.merge(mgf_info, cluster_info, how='left')

    def get_refined_dataset(self):

        """
        Return the refined dataset with updated cluster information.
        """

        return self.dataset

    def execute(self, search_ppm, rt):

        """
        Execute the cluster refinement process with given search_ppm and RT parameters.

        :param rt:
        :param search_ppm: m/z tolerance.
        """

        logger.info(f"Executing cluster refinement with search_ppm={search_ppm}, RT={rt}...")

        try:
            temp = self.dataset['Mass/Charge (m/z)'] * 1e-6 * search_ppm
            ppm_median = np.median(temp)

            self._refine_clusters_with_charge()  # Refine clusters based on charge
            self._ml_cluster_refinement('Cluster', ppm_median, rt, 1)  # Apply ML-based refinement

            logger.info("Cluster refinement executed successfully.")

        except Exception as e:
            logger.error(f"Error during execution: {e}")
            raise

    def _refine_clusters_with_charge(self):

        """
        Refine the dataset by redefining cluster indices based on charge.
        """

        logger.info("Refining clusters with charge...")
        try:
            # Combine cluster number with charge
            combined_clusters = self._combine_cluster_and_charge()

            # Map combined cluster numbers to integers
            self.dataset['Cluster'] = self._convert_to_integer_id(combined_clusters)

            logger.info("Clusters refined with charge successfully.")
        except Exception as e:
            logger.error(f"Error refining clusters with charge: {e}")
            raise

    def _combine_cluster_and_charge(self):

        """
        Combine the cluster number and charge to create unique cluster IDs.
        """

        return [f"{cluster}_{int(charge)}"
                for cluster, charge in self.dataset[['Cluster', 'Charge']].to_numpy()]

    @staticmethod
    def _convert_to_integer_id(combined_clusters):

        """
        Map the combined cluster IDs to unique integer values.
        """

        # Create a dictionary to map unique cluster IDs to integers
        unique_clusters = {cluster_id: idx for idx, cluster_id in enumerate(sorted(set(combined_clusters)))}

        # Map each cluster ID to its corresponding integer
        return [unique_clusters[cluster_id] for cluster_id in combined_clusters]

    def _ml_cluster_refinement(self, initial_cluster_column, ppm_eps, rt_ppm, min_samples):

        """
        Perform machine learning-based refinement of clusters.
        """

        logger.info("Starting ML-based cluster refinement...")
        try:
            top_1 = top1(self.dataset)
            top_1 = self._refine_clusters_with_ppm(top_1, initial_cluster_column, ppm_eps, min_samples)
            top_1 = self._refine_clusters_with_rt(top_1, 'New Cluster', rt_ppm, min_samples)
            top_1 = self._calculate_cluster_sizes(top_1, 'New Cluster')
            self._merge_cluster_info(top_1)

            logger.info("ML-based cluster refinement completed successfully.")

        except Exception as e:
            logger.error(f"Error during ML-based cluster refinement: {e}")
            raise

    @staticmethod
    def _refine_clusters_with_ppm(top_1, col, eps, min_samples):

        """
        Refine clusters based on Mass/Charge (m/z) and ppm.

        :param top_1:
        :param col: Column name to be used for refinement.
        :param eps: DBSCAN parameter for the maximum distance between samples.
        :param min_samples: DBSCAN parameter for the minimum number of samples in a cluster.
        :return: DataFrame with new refined cluster labels.
        """

        logger.info(f"Refining clusters with ppm. eps={eps}...")

        try:
            feature = top_1[[col, 'Mass/Charge (m/z)']]
            model = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1)
            top_1['New Cluster'] = model.fit_predict(feature)

            logger.info("Clusters refined with ppm successfully.")

            return top_1
        except Exception as e:
            logger.error(f"Error refining clusters with ppm: {e}")
            raise

    @staticmethod
    def _refine_clusters_with_rt(top_1, col, eps, min_samples):

        """
        Refine clusters based on retention time (RT).

        :param top_1:
        :param col: Column name to be used for refinement.
        :param eps: DBSCAN parameter for the maximum distance between samples.
        :param min_samples: DBSCAN parameter for the minimum number of samples in a cluster.
        :return: DataFrame with new refined cluster labels.
        """

        logger.info(f"Refining clusters with RT. eps={eps}...")

        try:
            feature = top_1[[col, 'Retention Time (min)']].copy()
            feature[col] = feature[[col]] * (eps + 1)  # Adjusting column by scaling with eps
            model = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1)
            top_1['New Cluster'] = model.fit_predict(feature)

            logger.info("Clusters refined with RT successfully.")

            return top_1
        except Exception as e:
            logger.error(f"Error refining clusters with RT: {e}")
            raise

    @staticmethod
    def _calculate_cluster_sizes(top_1, col):

        """
        Calculate the size of each cluster and add it as a new column.
        """

        logger.info("Calculating cluster sizes...")

        clu_list = top_1[col]
        clu_count = clu_list.value_counts().to_dict()
        top_1['New Count'] = top_1[col].map(clu_count)

        logger.info("Cluster sizes calculated successfully.")

        return top_1

    def _merge_cluster_info(self, top_1):

        """
        Merge the refined cluster information into the dataset.
        """

        self.dataset = pd.merge(self.dataset, top_1[['Source File', 'Scan number', 'New Cluster', 'New Count']],
                                on=['Source File', 'Scan number'], how='left')


class NewCandidatePeptideGenerator:
    def __init__(self, config, dataset):
        self.config = config
        self.dataset = dataset

    def select_top_n_candidates(self, top_n: int):

        """
        Selects the top n peptides based on scores for each cluster.

        Args:
            top_n (integer): Number of top peptides to select per cluster.

        Returns:
            defaultdict: {0: [((ABCD, 6.2345), 0),((ACBD, 5.64321), 1)]}
                        A dictionary where each cluster ID maps to a list of
                        tuples containing (peptide, score) and rank.
        """

        logging.info(f"Selecting top {top_n} candidates for each cluster...")

        self._remove_missing_values()
        score_dict = self._accumulate_scores()

        # Initializes a defaultdict to avoid manual key initialization
        top_n_peptides = defaultdict(list)

        for cluster_id in score_dict:
            # Sort peptides by score in descending order and take the top n
            sorted_peptides = sorted(score_dict[cluster_id].items(), key=lambda x: x[1], reverse=True)[:top_n]

            # Assign a rank starting from 0
            for rank, peptide_and_score in enumerate(sorted_peptides):
                top_n_peptides[cluster_id].append((peptide_and_score, rank))

        logging.info(f"Top {top_n} candidate selection completed successfully.")

        return top_n_peptides

    def _remove_missing_values(self):

        """
        Removes rows with missing values in 'Peptide' or 'Score' columns.

        Returns:
            pd.DataFrame: A DataFrame with missing values removed.
        """

        return self.dataset.dropna(subset=['Peptide', 'Score'])

    def _accumulate_scores(self):

        """
        Accumulates scores for each peptide within each cluster.
        """

        score_dict = {}

        for cluster_id, peptide, score in self.dataset[['New Cluster', 'Peptide', 'Score']].to_numpy():

            if cluster_id not in score_dict:
                score_dict[cluster_id] = {}

            if peptide not in score_dict[cluster_id]:
                score_dict[cluster_id][peptide] = 0
            score_dict[cluster_id][peptide] += score * 0.01

        # Previous version (commented out): Scores were stored in lists and summed later
        #     if peptide not in score_dict[cluster_id]:
        #         score_dict[cluster_id][peptide] = []
        #     score_dict[cluster_id][peptide].append(score * 0.01)
        #
        # for cluster_id in score_dict:
        #     for peptide in score_dict[cluster_id]:
        #         score_dict[cluster_id][peptide] = sum(score_dict[cluster_id][peptide])

        return score_dict

    def organize_new_candidates_info(self, top_peptides):

        """
        Organizes new candidate information by cluster, including peptides, sequences, ranks, and scores.

        Args:
            top_peptides (defaultdict): Top peptides with scores and ranks for each cluster.
       """

        logging.info("Organizing new candidate information by cluster...")

        peptides_by_cluster, sequences_by_cluster, ranks_by_cluster, scores_by_cluster = defaultdict(list), \
                                                                                         defaultdict(list), defaultdict(
            list), defaultdict(list)

        for peptide, sequence, cluster_id in self.dataset[['Peptide', 'Sequence', 'New Cluster']].to_numpy():

            for ranked_results in top_peptides[cluster_id]:
                if peptide == ranked_results[0][0]:
                    if peptide not in peptides_by_cluster[cluster_id]:
                        peptides_by_cluster[cluster_id].append(peptide)
                        sequences_by_cluster[cluster_id].append(sequence)
                        ranks_by_cluster[cluster_id].append(ranked_results[1])
                        scores_by_cluster[cluster_id].append(ranked_results[0][1])

        logging.info("Successfully organized new candidate information.")

        return peptides_by_cluster, sequences_by_cluster, ranks_by_cluster, scores_by_cluster

    @staticmethod
    def generate_new_candidates_dataframe(top_1_dataset, peptides_by_cluster, sequences_by_cluster, ranks_by_cluster,
                                          scores_by_cluster):

        """
        Generates a DataFrame of new candidate data.

        Args:
            top_1_dataset (pd.DataFrame): DataFrame containing the top 1 peptide data.

        Returns:
            pd.DataFrame: A DataFrame containing the new candidate peptides with their source file,
                          scan number, peptide, sequence, score, and rank.
        """

        logging.info("Generating new candidates DataFrame...")

        candidate_data = [
            (file_name, scan_number, peptide, sequence, score, rank)

            for file_name, scan_number, cluster_id in
            top_1_dataset[['Source File', 'Scan number', 'New Cluster']].to_numpy()

            for peptide, sequence, rank, score in
            zip(peptides_by_cluster[cluster_id], sequences_by_cluster[cluster_id],
                ranks_by_cluster[cluster_id], scores_by_cluster[cluster_id])]

        logging.info(f"New candidates DataFrame generated successfully. Total rows: {len(candidate_data)}")

        return pd.DataFrame(candidate_data,
                            columns=['Source File', 'Scan number', 'Peptide', 'Sequence', 'Score', 'Rank'])


class MGFNoiseRemover:
    def __init__(self, mgf_directory, save_directory):

        """
        Initialize the MGFNoiseRemover with the directories for MGF files and the denoised MGF save location.

        :param mgf_directory: Path to the directory containing the MGF files.
        :param save_directory: Path to the directory where the cleaned MGF files will be saved.
        """

        self.mgf_directory = mgf_directory
        self.save_directory = save_directory

    def generate_noise_removed_mgf(self):

        """
        Iterate through all MGF files in the input directory,
        process each spectrum to remove noise peaks, and save the denoised MGF files in the specified directory.
        """

        mgf_list = os.listdir(self.mgf_directory)

        for file_name in tqdm(mgf_list):
            self._process_mgf_file(file_name)

    def _process_mgf_file(self, filename):

        """
        Process a single MGF file: read its content, remove noise from spectra,
        and save the denoised data into a new file (with _remove added to the filename).

        :param filename: Name of the MGF file to process.
        """

        input_path = os.path.join(self.mgf_directory, filename)
        output_path = os.path.join(self.save_directory, f"{filename[:-4]}_remove.mgf")

        # Read input file content
        with open(input_path, 'r') as infile:
            data = infile.readlines()

        # Process the spectra and write denoised data to output file
        with open(output_path, 'w') as outfile:
            self._process_spectra(data, outfile)

    def _process_spectra(self, data, outfile):

        """
        Process the spectra data: remove noise, and write the denoised data to output.

        :param data: List of lines from the input MGF file.
        :param outfile: Output file object to write the denoised data.
        """

        peak_list = []
        window_min, window_max = 0, 100  # Initial window

        process_peaks = self._process_peaks
        adjust_window = self._adjust_window

        for line in data:
            # Start of a spectrum
            if line.startswith('BEGIN IONS'):
                # Reset the peak list for a new spectrum
                peak_list = []
                # Reset the window range
                window_min, window_max = 0, 100
                # Write the line to output file
                outfile.write(line)
            # End of a spectrum
            elif line.startswith('END IONS'):
                # Process the peaks
                process_peaks(peak_list, outfile)
                # Write the line to output file
                outfile.write(line)
            else:
                # Line contains a peak (mass-to-charge ratio and intensity)
                if line.strip() and line.split()[0].replace('.', '', 1).isdigit():
                    # Parse the mass-to-charge ratio (m/z)
                    mz = float(line.split()[0])
                    if mz <= window_max:
                        # Add the peak to the current list
                        peak_list.append(line)
                    else:
                        # Adjust the window range to include the new peak
                        window_min, window_max = adjust_window(mz, window_min, window_max)

                        # If there are no peaks in the current list, initialize it with the new peak
                        if len(peak_list) == 0:
                            peak_list = [line]
                            continue

                        # Process and write the current peaks, then reset for the new peak
                        process_peaks(peak_list, outfile)
                        peak_list = [line]

                # Non-peak lines (metadata) are directly written to the output
                else:
                    outfile.write(line)

    def _process_peaks(self, peak_list, outfile):

        """
        Process the peaks in the current spectrum: remove noise from the peaks within the m/z window range and write the denoised peaks.

        :param peak_list: List of peaks within the m/z window range in the current spectrum.
        :param outfile: Output file object to write the denoised peaks.
        """

        # Remove noise if there are more than 10 peaks
        if len(peak_list) > 10:
            peak_list = self._remove_noise_peak_keep_order(peak_list)

        # Write all peaks to the output file
        for peak in peak_list:
            outfile.write(peak)

    @staticmethod
    def _adjust_window(mz, window_min, window_max):

        """
        Adjust the m/z window range to accommodate a new peak.

        :param mz: Mass-to-charge ratio of the new peak.
        :param window_min: Current minimum of the window.
        :param window_max: Current maximum of the window.
        :return: Updated window_min and window_max.
        """

        # Increment the window by 100 until the m/z fits in the range
        while mz > window_max:
            window_min += 100
            window_max += 100

        return window_min, window_max

    @staticmethod
    def _remove_noise_peak_keep_order(peak_list):

        """
        Remove the least intense peaks one by one, ensuring only the top 10 peaks remain while maintaining the order of peaks.

        :param peak_list: List of peaks to process.
        :return: Denoised list of peaks with exactly 10 peaks.
        """

        intensity_list = [float(line.strip().split()[1]) for line in peak_list]

        # Keep removing the lowest intensity peak until only 10 peaks remain
        while len(peak_list) > 10:
            # Find the index of the lowest intensity
            min_index = intensity_list.index(min(intensity_list))

            # Remove the corresponding peak
            peak_list.pop(min_index)

            # Remove the intensity value
            intensity_list.pop(min_index)

        return peak_list


class XCorrMGFGenerator:
    def __init__(self, mgf_directory, save_directory):

        """
        Initialize the XcorrMGFGenerator.

        :param mgf_directory: Directory containing the original MGF files.
        :param save_directory: Directory to save the generated MGF files, modified for XCorr calculations.
        """

        self.mgf_directory = mgf_directory
        self.save_directory = save_directory

    def generate_cometx_mgf(self, dataset):

        """
        Generate new MGF files for XCorr calculation by adding peptide sequence information.
        """

        mgf_files = sorted(os.listdir(self.mgf_directory))  # List of MGF files sorted alphabetically
        current_file_index, current_index = 0, 0  # Track current MGF file and line index
        current_mgf_content, output_file = None, None  # Store MGF file content and output file reference
        mgf_scan_number = -1  # Current scan number being processed
        current_spectrum = []  # Store lines of the current spectrum

        for idx, (source_file, scan_number, peptide) in enumerate(
                tqdm(dataset[['Source File', 'Scan number', 'Peptide']].to_numpy(),
                     desc="Generating MGF files for XCorr calculation")):

            # Load the first MGF file at the start of processing
            if idx == 0:
                current_file_name = mgf_files[current_file_index]
                current_mgf_content = open(f"{self.mgf_directory}/{current_file_name}").readlines()
                output_file = open(os.path.join(self.save_directory, f"{current_file_name[:-4]}_xcorr.mgf"), 'w')
                current_index = 0

            if scan_number == mgf_scan_number:
                for data in current_spectrum:
                    if data.startswith('TITLE='):
                        formatted_peptide = peptide.replace('m', 'M+15.9949').replace('C', 'C+57.021464')
                        output_file.write(f'SEQ={formatted_peptide}\n')
                    output_file.write(data)
            else:
                # Process a new spectrum
                current_spectrum = []
                while True:
                    if current_index >= len(current_mgf_content):
                        current_file_index += 1
                        output_file.close()
                        if current_file_index >= len(mgf_files):  # If no more files, stop processing
                            logging.info("Reached the end of the files. No more data.")
                            break
                        # Load the next file
                        current_file_name = mgf_files[current_file_index]
                        current_mgf_content = open(f"{self.mgf_directory}/{current_file_name}").readlines()
                        output_file = open(os.path.join(self.save_directory, f"{current_file_name[:-4]}_xcorr.mgf"), 'w')
                        current_index = 0

                    current_line = current_mgf_content[current_index]

                    if current_line.startswith('TITLE='):
                        current_spectrum.append(current_line)
                        mgf_scan_number = int(current_line.split('.')[1])
                    elif current_line.startswith('END IONS'):
                        current_spectrum.append(current_line)
                        if mgf_scan_number == scan_number:
                            for data in current_spectrum:
                                if data.startswith('TITLE='):
                                    formatted_peptide = peptide.replace('m', 'M+15.9949').replace('C', 'C+57.021464')
                                    output_file.write(f'SEQ={formatted_peptide}\n')
                                output_file.write(data)
                            current_index += 1
                            break
                        else:
                            current_spectrum = []
                    else:
                        current_spectrum.append(current_line)
                    current_index += 1

        # Close the last output file
        output_file.close()
