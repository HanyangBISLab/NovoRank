import os
import re
import logging
import pandas as pd

from tqdm import tqdm
from src.utils.utils import *

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class FileHandler:
    def __init__(self, file_path: str):
        self.data = None
        self.file_path = file_path

    def load_csv(self) -> pd.DataFrame:

        """
        Load data from a CSV file in chunks.
        """

        if not self.file_path.endswith('.csv'):
            raise ValueError("Input file must be a .csv file.")

        logger.info(f"Loading CSV file: {self.file_path}")
        chunk_data = pd.read_csv(self.file_path, chunksize=100000)
        self.data = pd.concat(chunk_data)

        if 'Scan number' in self.data.columns:
            self.data = self.data.astype({'Scan number': 'int'})

        logger.info("CSV file loaded successfully.")
        return self.data


class SequencePreprocessor:
    def __init__(self):

        """
        ord(num): num is an uppercase letter if its value is between 65 and 91.
        """

        self.pattern = re.compile(r'[^0-9()+-.]')

    def process_denovo(self, de_novo) -> pd.DataFrame:

        """
        Process de novo sequences by cleaning and normalizing.
        """

        logger.info("Starting de novo sequence processing...")
        de_novo = de_novo[de_novo['Peptide'].notnull()]
        de_novo['Sequence'] = de_novo['Peptide'].apply(self._strip_sequence)

        logger.info("De novo sequence processing completed.")
        return de_novo.sort_values(by=['Source File', 'Scan number']).reset_index(drop=True)

    def process_database(self, db) -> pd.DataFrame:

        """
        Process database sequences by cleaning and normalizing.
        """

        logger.info("Starting database sequence processing...")
        db = db[db['GT'].notnull()]
        db['GT'] = db['GT'].apply(lambda seq: seq.replace('I', 'L'))

        logger.info("Database sequence processing completed.")
        return db.sort_values(by=['Source File', 'Scan number']).reset_index(drop=True)

    def _strip_sequence(self, sequence):

        """
        Normalize sequences by:
            - Converting isoleucine (I) to leucine (L).
            - Replacing oxidation-modified methionine (m) with methionine (M).
        """

        cleaned = ''.join(self.pattern.findall(sequence)).replace('I', 'L')
        cleaned = cleaned.replace('m', 'M')
        return cleaned


class ClusteringResultProcessor:
    def __init__(self, cluster_dir, mgf_dir):
        self.cluster_dir = cluster_dir
        self.cluster_list = os.listdir(cluster_dir)
        self.index_dic = self.create_file_index(os.listdir(mgf_dir))

    @staticmethod
    def create_file_index(mgf_list) -> dict:

        """
        Create an index dictionary from the MGF file list.
        """

        return {idx: mgf for idx, mgf in enumerate(mgf_list)}

    def process_clusters(self) -> pd.DataFrame:

        """
        Process the cluster files by reading them and organizing the cluster data.
        """

        logger.info("Starting cluster processing...")
        cluster = []

        for f in self.cluster_list:
            logger.info(f"Processing cluster file: {f}")
            with open(os.path.join(self.cluster_dir, f), 'r') as file:
                lines = file.readlines()

            clu_cnt, clu_num = 0, ""

            for line in lines:
                if not line.strip():  # Skip empty lines
                    continue
                parts = line.split('\t')

                if '.' in parts[0]:  # Line contains cluster information
                    clu_num = parts[0].split('.')[1]
                    clu_cnt = int(parts[1])
                else:  # Line contains spectrum information
                    idx = int(parts[1])
                    scan = int(parts[2])
                    fn = self.index_dic[idx]
                    cluster.append((fn, scan, clu_cnt, clu_num))

        logger.info("Cluster processing completed.")
        return pd.DataFrame(cluster, columns=['Source File', 'Scan number', 'Count', 'Cluster']).sort_values(['Source File', 'Scan number']).reset_index(drop=True)


class MGFProcessor:
    def __init__(self, mgf_path):
        self.path = mgf_path
        if not os.path.exists(mgf_path):
            raise FileNotFoundError(f"The specified path does not exist: {mgf_path}")

    def extract_spectrum_info(self) -> pd.DataFrame:

        """
        Extracts spectrum information from MGF files and returns a DataFrame.
        """

        mgf_files = [f for f in os.listdir(self.path) if f.endswith('.mgf')]
        if not mgf_files:
            logging.warning("No MGF files found in the specified path.")
            return pd.DataFrame()

        dataset = []
        logging.info("Starting spectrum information extraction...")

        for file_name in tqdm(mgf_files, desc="Processing MGF files"):
            with open(os.path.join(self.path, file_name), 'r') as file:
                scan, charge, mass, rt = None, None, None, None

                for line in file:
                    line = line.strip()
                    if line == 'END IONS':
                        if all(v is not None for v in [scan, charge, mass, rt]):
                            dataset.append((file_name, scan, charge, mass, rt))
                        else:
                            logging.warning(f"Incomplete spectrum data in file {file_name}")
                    elif line.startswith('TITLE='):
                        try:
                            scan = int(line.split('.')[1])
                        except (IndexError, ValueError):
                            logging.error(f"Failed to parse scan number from line: {line}")
                    elif line.startswith('CHARGE='):
                        try:
                            charge = int(line.split('=')[1].rstrip('+'))
                        except (IndexError, ValueError):
                            logging.error(f"Failed to parse charge from line: {line}")
                    elif line.startswith('RTINSECONDS='):
                        try:
                            rt = float(line.split('=')[1])
                        except (IndexError, ValueError):
                            logging.error(f"Failed to parse retention time from line: {line}")
                    elif line.startswith('PEPMASS='):
                        try:
                            mass = float(line.split('=')[1].split(' ')[0])
                        except (IndexError, ValueError):
                            logging.error(f"Failed to parse mass from line: {line}")

        columns = ['Source File', 'Scan number', 'Charge', 'Mass/Charge (m/z)', 'Retention Time (s)']
        mgf_info = pd.DataFrame(dataset, columns=columns)

        if not mgf_info.empty:
            mgf_info['Retention Time (min)'] = mgf_info['Retention Time (s)'] / 60
            mgf_info = mgf_info.drop(columns=['Retention Time (s)'])
            logging.info(f"Extracted spectrum information for {len(mgf_info)} spectra.")
        else:
            logging.warning("No valid information were extracted.")

        return mgf_info


class DataLoader:

    """
    Handles the loading and parsing of data for the training or testing process.

    Parameters:
        - config (dict): Configuration parameters.
    """

    def __init__(self, config):
        self.config = config

    def load_and_parse_data(self, cluster_info, mgf_info, ppm, rt) -> tuple[pd.DataFrame, pd.DataFrame]:

        """
        Returns:
            - de_novo (pd.DataFrame): Processed de_novo data.
            - clu (pd.DataFrame): Clustering result.
        """

        de_novo_path = os.path.join(self.config['path']['search_results']['de_novo']['path'],
                                    self.config['path']['search_results']['de_novo']['filename'])

        # Create instances of FileHandler and SequencePreprocessor for loading and preprocessing data
        file_handler = FileHandler(de_novo_path)
        preprocess = SequencePreprocessor()

        # Load and preprocess de_novo data
        de_novo = file_handler.load_csv()
        de_novo = preprocess.process_denovo(de_novo)
        refinement = ClusterRefinement(cluster_info, mgf_info)

        refinement.execute(ppm, rt)
        clu = refinement.get_refined_dataset()

        return de_novo, clu

    def load_db_data(self) -> pd.DataFrame:

        """
        Loads and preprocesses the database data for training.

        Returns:
            - db (pd.DataFrame): Processed database data.
        """

        db_path = os.path.join(self.config['path']['search_results']['db']['path'],
                               self.config['path']['search_results']['db']['filename'])

        # Create instances of FileHandler and SequencePreprocessor for loading and preprocessing data
        file_handler = FileHandler(db_path)
        preprocess = SequencePreprocessor()

        # Load and preprocess database data
        db = file_handler.load_csv()
        db = preprocess.process_database(db)

        return db


class CrossCorrelationResultProcessor:
    def __init__(self, xcorr_dir):

        """
        :param xcorr_dir: Path to the directory containing .tsv files with XCorr calculation results.
        """

        self.path_dir = xcorr_dir

    @staticmethod
    def _parse_line(line):

        """
        Parse a single line from the .tsv file.
        """

        line = line.strip()
        fn = line.split('.')[0] + '.mgf'
        scan = line.split('.')[1]
        charge = line.split('\t')[1].split(' ')[0]
        peptide = line.split('\t')[4].replace('M+15.9949', 'm').replace('C+57.021464', 'C')
        xcorr = line.split('\t')[5]

        return fn, scan, int(charge), peptide, float(xcorr)

    def extract_cross_correlation_info(self):

        """
        Process all the XCorr calculation result files in .tsv format and compile the data into a DataFrame.
        """

        file_ext = '.tsv'
        file_list = [file for file in os.listdir(self.path_dir) if file.endswith(file_ext)]

        data = []
        for file_name in tqdm(file_list, desc="Processing XCorr calculation result files"):
            file_path = os.path.join(self.path_dir, file_name)

            with open(file_path, 'r') as file:
                lines = file.readlines()

                for idx, line in enumerate(lines):
                    if idx == 0:  # Skip header
                        continue
                    parsed_data = self._parse_line(line)
                    data.append(parsed_data)

        columns = ['Source File', 'Scan number', 'Charge', 'Peptide', 'XCorr']
        df = pd.DataFrame(data, columns=columns).astype({'Scan number': 'int'})

        return df
