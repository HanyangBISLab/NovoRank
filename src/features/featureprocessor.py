import os
import logging
import numpy as np
import pandas as pd

from numba import jit
from tqdm import tqdm
from deeplc import DeepLC

# os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
# os.environ['CUDA_VISIBLE_DEVICES'] = '1'

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


amino_acid = {'G': 57.02146372376, 'A': 71.03711378804, 'S': 87.03202841014, 'P': 97.05276385232, 'V': 99.0684139166,
              'T': 101.04767847442, 'C': 160.03064447804, 'L': 113.08406398088, 'I': 113.08406398088, 'N': 114.04292744752,
              'D': 115.02694303224, 'Q': 128.0585775118, 'K': 128.09496301826, 'E': 129.04259309652, 'M': 131.0404846066,
              'H': 137.0589118628, 'F': 147.0684139166, 'R': 156.10111102874, 'Y': 163.0633285387, 'W': 186.07931295398,
              'm': 147.0353846066}

h20 = 18.01528
nh3 = 17.02655
proton = 1.00727647
neutron = 1.0086710869

fragment_tolerance = 0.025


class FeatureTransformer:
    def __init__(self, dataset):
        self.dataset = dataset

    def generate_features(self):

        """
        Generate new features for the dataset by applying log scaling and calculating delta scores.
        The results are added as new columns to the dataset.
        """

        # Apply log scaling to the dataset
        self._apply_log_scaling()

        # Calculate delta scores and add them as a new column 'Delta Score'
        self.dataset['Delta Score'] = self._calculate_delta_scores()

        # Reset the index
        return self.dataset.reset_index(drop=True)

    def _apply_log_scaling(self):

        """
        Apply logarithmic scaling to the 'Score' and 'New Count' columns to normalize the values.
        The transformed 'Score' is then divided by the transformed 'New Count' for scaling.
        """

        logging.info('Starting log scaling process...')

        # Apply log transformation to 'Score' and 'New Count'
        self.dataset['Score'] = np.log(self.dataset['Score'] + 1)
        self.dataset['New Count'] = np.log(self.dataset['New Count'] + 1)

        # Normalize 'Score' by dividing it by 'New Count'
        self.dataset['Score'] = self.dataset['Score'] / self.dataset['New Count']

        logging.info('Log scaling process completed successfully.')

    def _calculate_delta_scores(self):

        """
        Calculate the delta scores for each row.

        Returns:
            list: A list of delta scores corresponding to each row in the dataset.
        """

        logging.info('Starting delta score calculation process...')

        delta_scores = []  # List to store the calculated delta scores
        previous_key = None  # To store the previous key
        highest_score = 0  # To track the highest score

        for source_file, scan_number, score in self.dataset[['Source File', 'Scan number', 'Score']].values:
            current_key = (source_file, scan_number)

            if current_key == previous_key:
                # If still in the same group, calculate the delta score
                delta_scores.append(highest_score - score)
            else:
                # If a new group starts, initialize delta score to 0
                delta_scores.append(0)
                highest_score = score  # Reset highest score for the new group

            previous_key = current_key  # Update the current key
            highest_score = max(highest_score, score)  # Update the highest score

        logging.info('Delta score calculation process completed successfully.')

        return delta_scores


class featureExtractor:
    def __init__(self, config, elution_time):

        """
        Initialize the feature extractor with configuration settings.
        """

        self.denoised_mgf_path = config['path']['denoised_mgf_path']
        self.elution_time = elution_time

    def internal_fragment_ion_features(self, dataset):

        """
        Process the dataset to calculate internal fragment ion features.
        """

        # Separate data into those with missing sequences and those with sequences
        dataset_with_sequence = dataset[dataset['Sequence'].notnull()].copy()

        # Calculate internal fragment ions and sequence lengths
        dataset_with_sequence = self._calculate_internal_fragment_ion(dataset_with_sequence)
        dataset_with_sequence = self._calculate_sequence_length(dataset_with_sequence)

        # Normalize internal fragment ion counts by sequence length
        dataset_with_sequence['Normalized Internal Fragment Ions'] = dataset_with_sequence['Internal Fragment Ions Count'] / dataset_with_sequence['Sequence Length']
        dataset_with_sequence['Normalized Internal Fragment Ions'] = np.log(dataset_with_sequence['Normalized Internal Fragment Ions'] + 1)

        # Sort the processed data by specific columns and reset the index
        return dataset_with_sequence.sort_values(by=['Source File', 'Scan number', 'Rank']).reset_index(drop=True)

    def _calculate_internal_fragment_ion(self, dataset):

        """
       Calculate the count of internal fragment ions for each entry in the dataset
       using the provided denoised .mgf files.
       """

        process_ion_data = self._process_ion_data

        # Sort the list of .mgf files in the denoised path
        file_list = sorted(os.listdir(self.denoised_mgf_path))

        # Initialize variables
        current_file_index, line_index = 0, 0
        current_file_data, current_spec_key = None, None
        experimental_peak_list = []
        internal_fragment_ion = []

        # Iterate over each row of the dataset
        for source_file, scan_number, peptide, charge, mz_value in tqdm(
                dataset[['Source File', 'Scan number', 'Peptide', 'Charge', 'Mass/Charge (m/z)']].to_numpy(),
                desc="Finding internal fragment ions"):

            # Extract base name and scan number for spectrum identification
            file_base_name = source_file.split('.')[0]
            scan_number = int(scan_number)
            spec_key = (file_base_name, scan_number)

            # Load data if not already loaded, then initialize variables for the current file
            if current_file_data is None:
                current_file_name = file_list[current_file_index]
                current_file_data = open(f"{self.denoised_mgf_path}/{current_file_name}").readlines()
                line_index = 0

            # If processing the same spectrum, calculate and append the ion count
            if spec_key == current_spec_key:
                ion_count = process_ion_data(experimental_peak_list, charge, mz_value, peptide)
                internal_fragment_ion.append(ion_count)
                continue

            # Loop through the lines of the current .mgf file
            while True:
                # If the end of the file is reached, move to the next file
                if line_index >= len(current_file_data):
                    current_file_index += 1
                    if current_file_index >= len(file_list):  # No more files to process
                        logging.info("Reached the end of the files. No more data.")
                        break
                    current_file_name = file_list[current_file_index]
                    current_file_data = open(f"{self.denoised_mgf_path}/{current_file_name}").readlines()
                    line_index = 0

                # Process the current line in the .mgf file
                current_line = current_file_data[line_index].strip()  # Remove whitespace from both ends of the string

                if current_line.startswith('TITLE='):  # Identify spectrum titles

                    # Extract scan number and file name from the title
                    current_scan_number = int(current_line.split('.')[1])
                    current_file_name = current_line.split('.')[0].split('=')[1]
                    current_spec_key = (current_file_name, current_scan_number)

                    # Skip to the next file if current file is alphabetically less than the base name
                    if file_base_name > current_file_name:
                        current_file_index += 1
                        if current_file_index >= len(file_list):
                            logging.info("Reached the end of the files. No more data.")
                            break
                        current_file_name = file_list[current_file_index]
                        current_file_data = open(f"{self.denoised_mgf_path}/{current_file_name}").readlines()
                        line_index = 0
                        continue

                    # Initialize the experimental_peak_list
                    if spec_key == current_spec_key:
                        experimental_peak_list = []

                if spec_key == current_spec_key:
                    if current_line.startswith('END IONS'):  # Identify the end of spectrum data
                        # calculate and append the ion count
                        ion_count = process_ion_data(experimental_peak_list, charge, mz_value, peptide)
                        internal_fragment_ion.append(ion_count)
                        line_index += 1
                        break

                    # Collect peak data (mass/charge ratios)
                    if current_line.strip() and current_line.split()[0].replace('.', '', 1).isdigit():
                        mz = float(current_line.split()[0])
                        experimental_peak_list.append(mz)

                # Increment the line index
                line_index += 1

        # Add the calculated internal fragment ion counts to the dataset
        dataset['Internal Fragment Ions Count'] = internal_fragment_ion

        return dataset

    def _process_ion_data(self, experimental_peaks, charge, precursor_mass, peptide):

        """
        Processes internal fragment ion data for a given peptide.

        Parameters:
            experimental_peaks (list): List of m/z values from the spectrum.
            charge (integer): Charge state.
            precursor_mass (float): Precursor mass-to-charge ratio (m/z).
            peptide (string): Amino acid sequence of the peptide.

        Returns:
            int: Count of internal fragment ions.
        """

        # Calculate the b and y ion series masses for the given peptide sequence
        b_ion_list, y_ion_list = self._calculate_ion_series_masses(peptide)

        # Calculate the all possible theoretical ions
        theoretical_ion_list = self._calculate_residue_ions(b_ion_list, y_ion_list, int(charge), float(precursor_mass))

        # Remove peaks from the experimental data that match the calculated theoretical ions
        filtered_experimental_peaks = self._remove_matching_peaks(np.array(experimental_peaks), theoretical_ion_list)

        # Calculate the masses of all possible internal fragments of the peptide
        fragment_ion_mass_list = self._calculate_subpeptide_masses(peptide)

        # Consider charge effects on internal fragment ions
        internal_fragment_candidates = self._consider_charge(fragment_ion_mass_list, int(charge))

        # Count the final internal fragment ions present in the filtered experimental peaks
        ion_count = self._count_final_internal_fragments(filtered_experimental_peaks, internal_fragment_candidates)

        return ion_count

    @staticmethod
    def _calculate_sequence_length(dataset):

        """
        Calculate the length of the peptide sequence for each entry in the dataset.
        """

        # Compute the length of the peptide sequence for non-null entries
        dataset['Sequence Length'] = dataset['Sequence'].apply(lambda x: len(x))
        # apply(lambda x: len(x) if pd.notnull(x) else np.nan) -> Sequence: Nan, None 처리)

        return dataset

    @staticmethod
    @jit(nopython=False, cache=True)
    def _calculate_isotopes(masses):

        """
        Calculate isotopic shifts for a list of masses.

        Parameters:
            masses (numpy array): Array of base masses.

        Returns:
            numpy array: Array containing base masses along with +1 and +2 neutron shifts.
        """

        n = len(masses)
        isotopes = np.empty(n * 3)
        for i in range(n):
            base_mass = masses[i]

            # Include the base mass and isotopic shifts (+1 neutron, +2 neutrons)
            isotopes[i * 3] = base_mass
            isotopes[i * 3 + 1] = base_mass + neutron
            isotopes[i * 3 + 2] = base_mass + 2 * neutron

        return isotopes

    def _calculate_ion_series_masses(self, sequence):

        """
        Calculate b and y ion series masses for a peptide sequence.
        """

        # Map each amino acid in the sequence to its corresponding mass
        seq_masses = np.array([amino_acid[aa] for aa in sequence])
        # Reverse the sequence masses for y ion calculation
        rev_masses = seq_masses[::-1]

        # Compute cumulative sums to calculate b and y ion series masses
        # Calculate b-series masses
        b_masses = np.cumsum(seq_masses)
        # Calculate y-series masses (reverse of b-series)
        y_masses = np.cumsum(rev_masses)

        # Calculate isotopic masses
        b_isotopes = self._calculate_isotopes(b_masses)
        y_isotopes = self._calculate_isotopes(y_masses)

        # Sort and return results
        return np.sort(b_isotopes), np.sort(y_isotopes)

    @staticmethod
    @jit(nopython=False, cache=True)
    def _calculate_residue_ions(b, y, charge, precursor_mass):
        """
        Generate b and y ion series with isotopic charge states.

        Parameters:
        - b (numpy array): b-series masses.
        - y (numpy array): y-series masses.
        - charge (integer): maximum charge state to consider.
        - precursor_mass (float): precursor mass.

        Returns:
        - numpy array: unique mass values.
        """
        # Preallocate lists for b_ions and y_ions
        b_ions = []
        y_ions = []

        # Generate b-series ions
        for i in range(1, charge):
            b_ions.extend((b + (proton * i)) / i)
            b_ions.extend((b - h20 + (proton * i)) / i)
            b_ions.extend((b - nh3 + (proton * i)) / i)

        # Generate y-series ions
        for i in range(1, charge):
            y_ions.extend((y + h20 + (proton * i)) / i)
            y_ions.extend((y + (proton * i)) / i)  # h20 already included
            y_ions.extend((y - nh3 + h20 + (proton * i)) / i)

        # Include precursor mass
        y_ions.append(precursor_mass)

        # Combine and get unique results
        return np.unique(np.array(b_ions + y_ions))

    @staticmethod
    @jit(nopython=False, cache=True)
    def _remove_matching_peaks(experimental_peaks, theoretical_peaks):

        """
        Remove peaks from the experimental data that match theoretical peaks.

        Parameters:
        - experimental_peaks (numpy array): observed experimental peaks.
        - theoretical_peaks (numpy array): calculated theoretical peaks.

        Returns:
        - numpy array: filtered experimental peaks.
        """

        # Initialize a mask array with all values set to True (indicating that all peaks are valid initially)
        mask = np.ones(len(experimental_peaks), dtype=np.bool_)

        # Check if the peak is within the tolerance range of any theoretical peak
        for idx, val in enumerate(experimental_peaks):
            for i in theoretical_peaks:
                if i - fragment_tolerance <= val <= i + fragment_tolerance:
                    mask[idx] = False
                    break

        # Return the experimental peaks that are still marked as valid (True)
        return experimental_peaks[mask]

    @staticmethod
    def _calculate_subpeptide_masses(peptide):

        """
        Calculate masses for all theoretical internal fragment ions.
        """

        # Remove the first and last character of the peptide
        subpeptide_sequence = peptide[1:-1]

        # Generate all possible subpeptides of varying lengths
        subpeptide_set = {subpeptide_sequence[i:i + n]
                          for n in range(1, len(subpeptide_sequence) + 1)
                          for i in range(len(subpeptide_sequence) - n + 1)}

        # Calculate the mass of each subpeptide
        subpeptide_masses = {sum(amino_acid[aa] for aa in subpeptide)
                             for subpeptide in subpeptide_set}

        return np.array(sorted(subpeptide_masses))

    @staticmethod
    def _consider_charge(fragment_ion_mass_list, charge):

        """
        Calculate ion masses for multiple charge states.

        Parameters:
        - fragment_ion_mass_list (numpy array): list of base masses.
        - charge (integer): charge state.

        Returns:
        - numpy array: unique ion masses.
        """

        # compute fragment ions masses for all charge states dynamically
        ion_masses = [(fragment_ion_mass_list + (proton * i)) / i for i in range(1, charge)]

        # Convert to a numpy array and remove duplicates
        return np.unique(np.array(ion_masses))

    @staticmethod
    @jit(nopython=False, cache=True)
    def _count_final_internal_fragments(filtered_experiment_peaks, theoretical_internal_fragments):

        """
        Count the number of matching internal fragment peaks.

        Parameters:
        - filtered_experiment_peaks (numpy array): experimental peaks after filtering.
        - theoretical_internal_fragments (numpy array): theoretical internal fragment peaks.

        Returns:
        - integer: count of matching internal fragment peaks.
        """

        # Initialize a mask array to track which peaks match the theoretical fragments
        mask = np.zeros(len(filtered_experiment_peaks), dtype=np.bool_)

        # Check if the peak is within the tolerance range of any theoretical internal fragment
        for idx, peak in enumerate(filtered_experiment_peaks):
            for i in theoretical_internal_fragments:
                if i - fragment_tolerance <= peak <= i + fragment_tolerance:
                    mask[idx] = True
                    break

        # Count and return the number of matching peaks
        return np.count_nonzero(mask)

    def calculate_retention_time_difference_features(self, dataset):

        """
        Calculates the difference between the predicted and observed retention times.
        """

        # Normalize retention time (in minutes)
        dataset['Retention Time (min)'] = dataset['Retention Time (min)'] / (self.elution_time / 60)

        # Data with valid peptide sequences
        valid_sequence_data = dataset[dataset['Sequence'].notnull()].copy()

        # Select the top two peptides for each scan(=spectrum)
        valid_sequence_data = self._select_top_two_peptides(valid_sequence_data)

        # Prepare specific data for retention time prediction
        peptide_data, calibration_data = self._prepare_peptide_and_calibration_data(valid_sequence_data)

        # Predict retention times
        predicted_rt = self._predict_retention_times(peptide_data, calibration_data)

        # Add predicted retention times to the valid sequence data
        valid_sequence_data['Predicted_RT (min)'] = predicted_rt

        # Calculate the logarithmic absolute difference between observed and predicted retention times
        valid_sequence_data['Difference_RT (min)'] = np.log(abs(valid_sequence_data['Retention Time (min)'] - valid_sequence_data['Predicted_RT (min)']))

        # Sort the processed data by specific columns and reset the index
        return valid_sequence_data.sort_values(by=['Source File', 'Scan number', 'Rank']).reset_index(drop=True)

    @staticmethod
    def _select_top_two_peptides(valid_sequence_data):

        """
        Selects the top two peptides for each spectrum.
        """

        dataset = valid_sequence_data.sort_values(by=['Source File', 'Scan number', 'Rank'])

        # Filter out peptides with lengths greater than or equal to 60
        dataset = dataset[dataset['Peptide'].str.len() < 60].reset_index(drop=True)

        # For each spectrum, select the first (highest rank) and last (lowest rank) entries
        first_rank = dataset.drop_duplicates(subset=['Source File', 'Scan number'], keep='first')
        last_rank = dataset.drop_duplicates(subset=['Source File', 'Scan number'], keep='last')

        # Combine these entries into a single dataset
        return pd.concat([first_rank, last_rank]).sort_values(
            by=['Source File', 'Scan number', 'Rank']
        ).reset_index(drop=True)

    @staticmethod
    def _prepare_peptide_and_calibration_data(valid_sequence_data):

        """
        Prepares data for retention time prediction.
        """

        # Extract relevant columns
        dataset = valid_sequence_data[['Source File', 'Scan number', 'Sequence', 'Peptide', 'Retention Time (min)', 'Score']].copy()

        # Convert modifications to strings in a specific format
        modifications = []
        for peptide in dataset['Peptide']:
            mods = [
                f"{idx + 1}|Carbamidomethyl" if aa == "C" else f"{idx + 1}|Oxidation"
                for idx, aa in enumerate(peptide)
                if aa in ['C', 'm']
            ]
            modifications.append("|".join(mods))
        dataset['modifications'] = modifications

        # Rename columns
        dataset.rename(columns={'Retention Time (min)': 'tr', 'Sequence': 'seq'}, inplace=True)

        # Split the dataset
        cal_data = dataset[['Source File', 'Scan number', 'seq', 'modifications', 'tr', 'Score']].copy()
        cal_data['modifications'] = ""

        pep_data = dataset[['Source File', 'seq', 'modifications', 'tr']].copy()

        # Fill missing modification values with an empty string
        cal_data['modifications'] = cal_data['modifications'].fillna("")
        pep_data['modifications'] = pep_data['modifications'].fillna("")

        return pep_data, cal_data

    def _predict_retention_times(self, peptide_data, calibration_data):

        """
        Predicts retention times for peptides using DeepLC.
        """

        predicted_rt = []

        # Importing DeepLC
        dlc = DeepLC()

        for source_file in tqdm(sorted(peptide_data['Source File'].unique()), desc="Predicting retention times"):
            peptide_subset = peptide_data[peptide_data['Source File'] == source_file][['seq', 'modifications', 'tr']]
            calibration_subset = calibration_data[calibration_data['Source File'] == source_file]

            # Prepare calibration data
            calibration_subset = self._prepare_calibration_data(calibration_subset)

            # Perform the calibration process
            dlc.calibrate_preds(seq_df=calibration_subset)

            # Use DeepLC to predict retention times for the input peptides
            preds = dlc.make_preds(seq_df=peptide_subset)
            predicted_rt.extend(preds)

        return predicted_rt

    @staticmethod
    def _prepare_calibration_data(dataset):

        """
        Prepares calibration by filtering the data.
        """

        # Remove duplicate entries
        dataset = dataset.drop_duplicates(subset=['Source File', 'Scan number'], keep='first')
        dataset = dataset.drop_duplicates(subset=['seq', 'modifications'], keep='first')

        # Select the top 1000 entries based on score
        dataset = dataset.sort_values(by=['Score'], ascending=False).head(1000)

        return dataset[['seq', 'modifications', 'tr']]
