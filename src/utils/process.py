import shutil
import logging
import pandas as pd

from src.utils.utils import *
from src.loader.dataloader import *

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class MergeProcess:
    def __init__(self, de_novo, cluster_info, db=None):
        self.de_novo = de_novo
        self.cluster_info = cluster_info
        self.db = db

    def merge_data(self) -> pd.DataFrame:

        """
        Merges the de novo, clustering, and optionally database search results.

        Returns:
            pd.DataFrame: Final merged and cleaned dataset.
        """

        logger.info("Starting the merging process for de novo, clustering, and database search results...")
        if self.db is not None:
            # Outer merge between de novo and database search results
            dataset = pd.merge(self.de_novo, self.db, how='outer')

            # Left merge with clustering results
            dataset = pd.merge(dataset, self.cluster_info, how='left')

            # Retain only the top-ranked result
            # Count and log the number of reliable PSMs (GT)
            logger.info(f"Reliable PSM (GT): {top1(dataset)['GT'].notnull().sum()} scans.")

            # Remove rows with missing clustering information
            dataset = drop_missing_cluster_info(dataset)
        else:
            # Left merge with clustering results
            dataset = pd.merge(self.de_novo, self.cluster_info, how='left')

            # Remove rows with missing clustering information
            dataset = drop_missing_cluster_info(dataset)

        logger.info("Merging process completed successfully.")
        return dataset.reset_index(drop=True)


class TrainProcess:

    """
    Handles the end-to-end data preparation for the 'NovoRank' training process.

    Responsibilities:
        - Load and parse data from input datasets.
        - Merge these datasets into a single dataframe for training.

    Parameters:
        - config (dictionary): Configuration settings (e.g., loaded from config.yaml).
    """

    def __init__(self, config):
        self.data_loader = DataLoader(config)
        self.config = config

    def execute_data_processing(self, cluster_info, mgf_info, ppm, rt) -> pd.DataFrame:

        """
        Execute the data loading, merging, and preparation pipeline for training.

        Parameters:
            - cluster_info (pd.DataFrame): Clustering results as a dataframe.
            - mgf_info (pd.DataFrame): Spectrum information as a dataframe.
        """

        # Load and parse de novo and clustering data
        de_novo, cluster = self.data_loader.load_and_parse_data(cluster_info, mgf_info, ppm, rt)

        # Load and parse database data
        db = self.data_loader.load_db_data()

        # Merge datasets using the MergeProcess class
        merge_obj = MergeProcess(de_novo, cluster, db)
        merged_df = merge_obj.merge_data()

        return merged_df

    @staticmethod
    def execute_candidate_generation(dataset, dataset_top1):

        """
        This function merges the input dataset with `dataset_top1`, performs filtering and sorting to generate the final new dataset.
        """

        # Select specific columns from the dataset_top1 for merging
        dataset_top1 = dataset_top1[['Source File', 'Scan number', 'Charge', 'Mass/Charge (m/z)', 'Retention Time (min)', 'New Cluster', 'New Count', 'GT']]

        # Merge the new dataset with the dataset_top1 from original based on common columns
        dataset = pd.merge(dataset, dataset_top1, how='outer')

        # Apply the labeling function to add a label column based on peptide and GT match
        dataset = labeling(dataset)

        # Remove rows where 'Peptide' is missing
        dataset = dataset[dataset['Peptide'].notnull()]

        # Sort the dataset by 'Source File', 'Scan number', and 'Rank', and reset the index
        return dataset.sort_values(by=['Source File', 'Scan number', 'Rank']).reset_index(drop=True)

    def execute_xcorr_mgf_generation(self, dataset):

        """
        Execute the process of generating MGF files for XCorr calculation, filtering the dataset beforehand.
        """

        # Initialize the process with configuration settings
        process = XCorrMGFGenerateProcess(self.config)
        # Set up a clean output directory
        process.prepare_directory()

        # Filter rows with non-null 'GT' values
        filtered_dataset = dataset[dataset['GT'].notnull()].reset_index(drop=True)

        # Generate an MGF files containing the peptide sequence information
        process.execute(filtered_dataset)

        return filtered_dataset


class InferenceProcess:

    """
    Handles the end-to-end data preparation for the 'NovoRank' inference process.

    Responsibilities:
        - Load and parse data from input datasets.
        - Merge these datasets into a single dataframe for inference.

    Parameters:
        - config (dictionary): Configuration settings (e.g., loaded from config.yaml).
    """

    def __init__(self, config):
        self.data_loader = DataLoader(config)

    def execute_data_processing(self, cluster_info, mgf_info, ppm, rt) -> pd.DataFrame:

        """
        Executes the process of loading, parsing, and merging datasets for the test phase.

        Parameters:
            - cluster_info (pd.DataFrame): Clustering results as a dataframe.
            - mgf_info (pd.DataFrame): Spectrum information as a dataframe.
        """

        # Load and parse de novo and clustering data
        de_novo, cluster = self.data_loader.load_and_parse_data(cluster_info, mgf_info, ppm, rt)

        # Merge de novo and clustering results (no DB data in testing phase)
        merge_obj = MergeProcess(de_novo, cluster)
        merged_df = merge_obj.merge_data()

        return merged_df

    @staticmethod
    def execute_candidate_generation(dataset, dataset_top1):

        """
        This function merges the input dataset with `dataset_top1`, performs filtering and sorting to generate the final new dataset.
        """

        # Select specific columns from dataset_top1
        dataset_top1 = dataset_top1[['Source File', 'Scan number', 'Charge', 'Mass/Charge (m/z)', 'Retention Time (min)', 'New Cluster', 'New Count']]

        # Merge new dataset with the dataset_top1 from original using an outer join
        dataset = pd.merge(dataset, dataset_top1, how='outer')

        # Remove rows where 'Peptide' is missing
        dataset = dataset[dataset['Peptide'].notnull()]

        # Sort the dataset by 'Source File', 'Scan number', and 'Rank', and reset the index
        return dataset.sort_values(by=['Source File', 'Scan number', 'Rank']).reset_index(drop=True)

    def execute_xcorr_mgf_generation(self, dataset):

        """
        Execute the entire process of generating MGF files for XCorr calculation.
        """

        # Initialize the process with configuration settings
        process = XCorrMGFGenerateProcess(self.config)
        # Set up a clean output directory
        process.prepare_directory()

        # Generate an MGF files containing the peptide sequence information
        process.execute(dataset)

        return dataset


class GenerateNewCandidateProcess:

    """
    Handles the generation of new peptide candidates from a given dataset.

    Key Responsibilities:
        - Selects the top N peptides per cluster to generate new peptide candidates.
        - Collects and organizes peptide, sequence, rank, and score information.
        - Creates a dataset with the newly selected top N candidates.

    Parameters:
        - config (dictionary): Configuration settings (e.g., loaded from config.yaml).
        - dataset (pd.DataFrame): The dataset containing peptide and other associated information.
    """

    def __init__(self, config, dataset):
        self.config = config
        self.dataset = dataset

        # Initialize PeptideCandidateGenerator
        self.peptide_generator = NewCandidatePeptideGenerator(config, dataset)

    def execute(self, top_n: int) -> pd.DataFrame:

        """
        Executes the new peptide candidates generation process.

        Parameters:
            - top_n (integer): The number of top peptide candidates to select.
        """

        # Generate top N candidates
        n_candidates = self.peptide_generator.select_top_n_candidates(top_n)

        # Organize peptides, sequences, ranks, and scores for the top N candidates
        peptides, sequences, ranks, scores = self.peptide_generator.organize_new_candidates_info(n_candidates)

        # Prepare the dataset with top 1 candidates
        dataset_top1 = top1(self.dataset)

        # Create a dataset consisting of newly defined peptide candidates
        new_dataset = self.peptide_generator.generate_new_candidates_dataframe(dataset_top1, peptides, sequences, ranks, scores)

        return new_dataset, dataset_top1


class SpectrumNoiseRemoveProcess:
    def __init__(self, config):

        """
        Initialize the SpectrumNoiseRemoverProcess with configuration settings.
        """

        self.mgf_directory = config['path']['mgf_path']  # Path to the directory containing original MGF files
        self.save_directory = config['path']['denoised_mgf_path']  # Path to the directory for saving denoised MGF files
        self.instance = MGFNoiseRemover(self.mgf_directory, self.save_directory)  # Create an MGFNoiseRemover instance

    def prepare_directory(self):

        """
        Prepare the output directory for saving denoised MGF files.

        If the directory already exists, it is removed and recreated to ensure a clean state.
        """

        if os.path.exists(self.save_directory):  # Check if the output directory exists
            shutil.rmtree(self.save_directory)  # Remove the directory if it exists
        os.mkdir(self.save_directory)  # Create a fresh output directory

    def execute(self):

        """
        Execute the noise removal process for all original MGF files.

        Saves the denoised spectra to the specified output directory.
        """

        # denoised MGF files: Each file contains spectra where only the top 10 peaks are retained within every 100 Da window
        self.instance.generate_noise_removed_mgf()


class XCorrMGFGenerateProcess:
    def __init__(self, config):

        """
        Initialize the XcorrMGFGenerator with configuration settings.
        """

        self.mgf_directory = config['path']['mgf_path']  # Path to the directory containing original MGF files
        self.save_directory = config['path']['xcorr_mgf_path']  # Path to the directory for saving MGF files for Xcorr calculation
        self.instance = XCorrMGFGenerator(self.mgf_directory, self.save_directory)  # Create an XCorrMGFGenerator instance

    def prepare_directory(self):

        """
        Prepare the output directory for saving processed MGF files.

        If the directory already exists, it is removed and recreated to ensure a clean state.
        """

        if os.path.exists(self.save_directory):  # Check if the output directory exists
            shutil.rmtree(self.save_directory)  # Remove the directory if it exists
        os.mkdir(self.save_directory)  # Create a fresh output directory

    def execute(self, dataset):

        """
        Execute the MGF file generation process for XCorr calculation.

        Generates new MGF files containing the peptide sequence information, saving the output in the specified directory.
        """

        self.instance.generate_cometx_mgf(dataset)
