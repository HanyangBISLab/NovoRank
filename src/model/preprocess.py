import os
import logging
import numpy as np
import pandas as pd

from numba import jit
from sklearn.model_selection import GroupShuffleSplit

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class DataPreparation:

    """
    A class to prepare datasets for training and inference.
    Includes methods for splitting datasets into training, validation, and test sets.
    """

    def __init__(self):
        self.valid_df = None
        self.above_max_seq_df = None
        self.missing_xcorr_df = None

    def get_train_dataset(self, dataset: pd.DataFrame, xcorr_info: pd.DataFrame):

        """
        Prepare the dataset for training by retaining valid sequences.
        """

        # Merge the dataset with XCorr information
        merged_df = pd.merge(
            dataset, xcorr_info,
            on=['Source File', 'Scan number', 'Peptide', 'Charge'],
            how='outer'
        )
        # merged_df['XCorr'] = np.log(1 + merged_df['XCorr'])

        # Filter out rows where Sequence or XCorr values are missing
        filtered_df = merged_df[merged_df['Sequence'].notnull() & merged_df['XCorr'].notnull()]
        filtered_df = filtered_df.reset_index(drop=True)

        # Remove duplicates and calculate Delta XCorr
        processed_df = self._remove_duplicates(filtered_df)
        processed_df['Delta XCorr'] = processed_df['XCorr_x'] - processed_df['XCorr_y']

        # Retain only sequences with acceptable lengths
        self.valid_df = self._get_length_filtered_data(processed_df)

        # Use only data with a clear label (one label is 1, the other is 0)
        self.valid_df = self.valid_df[(self.valid_df['Label_x'] + self.valid_df['Label_y']) == 1]

        return self.valid_df

    def get_inference_dataset(self, dataset: pd.DataFrame, xcorr_info: pd.DataFrame):

        """
        Prepares the dataset for inference, including remaining datasets.
        """

        # Merge the dataset with XCorr information
        merged_df = pd.merge(
            dataset, xcorr_info,
            on=['Source File', 'Scan number', 'Peptide', 'Charge'],
            how='outer'
        )
        # merged_df['XCorr'] = np.log(1 + merged_df['XCorr'])

        # Filter out rows where Sequence or XCorr values are missing
        filtered_df = merged_df[merged_df['Sequence'].notnull() & merged_df['XCorr'].notnull()]
        filtered_df = filtered_df.reset_index(drop=True)

        # Remove duplicates and calculate Delta XCorr
        processed_df = self._remove_duplicates(filtered_df)
        processed_df['Delta XCorr'] = processed_df['XCorr_x'] - processed_df['XCorr_y']

        # Separate valid sequences and sequences exceeding maximum length
        self.valid_df, self.above_max_seq_df = self._get_length_filtered_data(processed_df, False)

        # Rows with XCorr values missing but sequences available
        self.missing_xcorr_df = merged_df[merged_df['Sequence'].notnull()
                                          & merged_df['XCorr'].isnull()].reset_index(drop=True)

        return self.valid_df, self.above_max_seq_df, self.missing_xcorr_df

    @staticmethod
    def _remove_duplicates(dataset: pd.DataFrame) -> pd.DataFrame:

        """
        Retains only the first and last rows for each key defined by 'Source File' and 'Scan number',
        removing all other rows.
        """

        # Sort dataset by source file, scan number, and rank
        sorted_df = dataset.sort_values(by=['Source File', 'Scan number', 'Rank'])

        # Retain the first occurrence of each group
        top_df = sorted_df.drop_duplicates(subset=['Source File', 'Scan number'], keep='first')
        # Retain the last occurrence of each group
        bottom_df = sorted_df.drop_duplicates(subset=['Source File', 'Scan number'], keep='last')

        # Merge the first and last rows of each group
        merged_df = pd.merge(
            top_df, bottom_df,
            on=['Source File', 'Scan number'],
            how='outer'
        ).reset_index(drop=True)

        return merged_df

    @staticmethod
    def _get_length_filtered_data(dataset: pd.DataFrame, train=True):

        """
        Filters the dataset based on the length of peptide sequences.
        For training, retains only sequences within acceptable length bounds.
        For inference, also separates sequences exceeding the length bounds.
        """

        # Calculate the lengths of the peptide sequences
        dataset['Length_x'] = dataset['Peptide_x'].apply(len)
        dataset['Length_y'] = dataset['Peptide_y'].apply(len)

        # Identify sequences within acceptable length bounds
        peptide_below_x = dataset['Length_x'] <= 40
        peptide_below_y = dataset['Length_y'] <= 40

        dataset_below = dataset[peptide_below_x & peptide_below_y].reset_index(drop=True)

        if train:
            return dataset_below
        else:
            # Identify sequences exceeding the maximum length bounds
            peptide_above_x = dataset['Length_x'] > 40
            peptide_above_y = dataset['Length_y'] > 40

            dataset_above = dataset[peptide_above_x | peptide_above_y].reset_index(drop=True)

            return dataset_below, dataset_above

    @staticmethod
    def train_val_split(dataset, val_size, test=False):

        """
        Splits the dataset into training and validation sets (and optionally test sets) using group-based splitting.
        """

        # Group by ground truth labels
        groups = dataset['GT_x']
        gss = GroupShuffleSplit(n_splits=1, test_size=val_size)

        # Get train and validation indices
        train_idx, val_idx = next(gss.split(dataset, groups=groups))

        train_data = dataset.iloc[train_idx].reset_index(drop=True)
        val_data = dataset.iloc[val_idx].reset_index(drop=True)

        # Log the size of the training set
        logger.info(f"Train_data_size: {len(train_data)}")

        # If test is True, further split the validation set into validation and test sets
        if test:
            # Further split the validation set into validation and test sets
            val_gss = GroupShuffleSplit(n_splits=1, test_size=0.5)
            _, test_idx = next(val_gss.split(val_data, groups=val_data['GT_x']))

            # Extract the data using the indices
            test_data = val_data.iloc[test_idx].reset_index(drop=True)
            val_data = val_data.drop(test_idx).reset_index(drop=True)

            # Log the size of validation and test sets
            logger.info(f"Validation_data_size: {len(val_data)}")
            logger.info(f"Test_data_size: {len(test_data)}")

            return train_data, val_data, test_data
        else:
            # Return train and validation sets without splitting into test set
            # Log the size of the validation set
            logger.info(f"Validation_data_size: {len(val_data)}")
            return train_data, val_data


class InputDataGenerator:
    def __init__(self, mgf_path, max_length=50000, mz_resolution=0.1):
        self.mgf_path = mgf_path
        self.max_length = max_length
        self.mz_resolution = mz_resolution

        # Mapping table for charge states and amino acids for sequence encoding
        self.charge_index = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5}
        self.amino_acid_index = {
            'A': 6, 'C': 7, 'D': 8, 'E': 9, 'F': 10, 'G': 11, 'H': 12, 'I': 13, 'K': 14, 'L': 15,
            'M': 16, 'm': 16, 'N': 17, 'P': 18, 'Q': 19, 'R': 20, 'S': 21, 'T': 22, 'V': 23, 'W': 24, 'Y': 25
        }

        # Modifications for specific residues
        self.modification_index = {'m': 26, 'C': 27}

    def generate_data(self, dataset):

        """
        Generates input data for model training, validation by processing spectrum files and metadata.
        """

        # Initialize variables for file processing
        current_file_index, line_index = 0, 0
        current_file_data, current_spec_key = None, None
        spectrum = None

        # Sort the list of .mgf files in the directory
        file_list = sorted(os.listdir(self.mgf_path))

        # Generate input data for deep learning
        for (
                source_file, scan_number, peptide_x, peptide_y, charge, cscore_x, cscore_y, delta_cscore,
                ifi_x, ifi_y, diff_rt_x, diff_rt_y, xcorr_x, xcorr_y, delta_xcorr, label
            ) in dataset[['Source File', 'Scan number', 'Peptide_x', 'Peptide_y',
                          'Charge_x', 'Score_x', 'Score_y', 'Delta Score_y',
                          'Normalized Internal Fragment Ions_x', 'Normalized Internal Fragment Ions_y',
                          'Difference_RT (min)_x', 'Difference_RT (min)_y',
                          'XCorr_x', 'XCorr_y', 'Delta XCorr', 'Label_x']].values:

            # Extract spectrum identification details
            file_base_name = source_file.split('.')[0]
            scan_number = int(scan_number)
            spec_key = (file_base_name, scan_number)

            # Load data if not already loaded
            if current_file_data is None:
                current_file_name = file_list[current_file_index]
                with open(f"{self.mgf_path}/{current_file_name}") as f:
                    current_file_data = f.readlines()
                line_index = 0

            while True:
                # If the current file has reached the end
                if line_index >= len(current_file_data):
                    current_file_index += 1

                    # No more files to process
                    if current_file_index >= len(file_list):
                        logging.info("Reached the end of the files. No more data.")
                        break
                    current_file_name = file_list[current_file_index]
                    with open(f"{self.mgf_path}/{current_file_name}") as f:
                        current_file_data = f.readlines()
                    line_index = 0

                # Read the current line
                current_line = current_file_data[line_index].strip()

                if current_line.startswith('TITLE='):
                    # Extract scan number and file name from the TITLE
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
                        with open(f"{self.mgf_path}/{current_file_name}") as f:
                            current_file_data = f.readlines()
                        line_index = 0
                        continue

                    # Initialize the spectrum if the current spectrum matches
                    if spec_key == current_spec_key:
                        spectrum = self._initialize_spectrum(self.max_length)

                if spec_key == current_spec_key:
                    if current_line.startswith('END IONS'):
                        # Normalize and encode the spectrum and metadata
                        spectrum = self._normalize_spectrum(spectrum)
                        sequence_x = self._encode_sequence(peptide_x, charge)
                        sequence_y = self._encode_sequence(peptide_y, charge)
                        psm_features_x = self._add_psm_features(cscore_x, ifi_x, diff_rt_x, xcorr_x)
                        psm_features_y = self._add_psm_features(cscore_y, ifi_y, diff_rt_y, xcorr_y)
                        delta_features = self._add_delta_features(delta_cscore, delta_xcorr)  # count
                        label_ = self._labeling(label)

                        # Yield the processed input data and label
                        yield {
                            'input_spectrum': spectrum,
                            'input_sequence_x': sequence_x,
                            'input_sequence_y': sequence_y,
                            'input_psm_features_x': psm_features_x,
                            'input_psm_features_y': psm_features_y,
                            'input_delta_features': delta_features
                        }, label_

                        line_index += 1
                        break

                    # Collect peak data if it's a valid spectrum line
                    if current_line.strip() and current_line.split()[0].replace('.', '', 1).isdigit():
                        mz = float(current_line.split()[0])
                        intensity = float(current_line.split()[1])
                        loc = int(mz // self.mz_resolution)
                        if loc <= self.max_length:
                            spectrum[loc] += intensity

                # Move to the next line
                line_index += 1

    def generate_inference_data(self, dataset, input_names):

        """
        Generates input data for inference by processing spectrum files and metadata.
        """

        # Unpack input names for clarity and convenience
        input_spectrum = input_names[0]
        input_sequence_x = input_names[1]
        input_sequence_y = input_names[2]
        input_psm_features_x = input_names[3]
        input_psm_features_y = input_names[4]
        input_delta_features = input_names[5]

        # Initialize variables for file processing
        current_file_index, line_index = 0, 0
        current_file_data, current_spec_key = None, None
        spectrum = None

        # Sort the list of .mgf files in the directory
        file_list = sorted(os.listdir(self.mgf_path))

        # Generate input data for deep learning
        for (
                source_file, scan_number, peptide_x, peptide_y, charge, cscore_x, cscore_y, delta_cscore,
                ifi_x, ifi_y, diff_rt_x, diff_rt_y, xcorr_x, xcorr_y, delta_xcorr
            ) in dataset[['Source File', 'Scan number', 'Peptide_x', 'Peptide_y',
                          'Charge_x', 'Score_x', 'Score_y', 'Delta Score_y',
                          'Normalized Internal Fragment Ions_x', 'Normalized Internal Fragment Ions_y',
                          'Difference_RT (min)_x', 'Difference_RT (min)_y',
                          'XCorr_x', 'XCorr_y', 'Delta XCorr']].values:

            # Extract spectrum identification details
            file_base_name = source_file.split('.')[0]
            scan_number = int(scan_number)
            spec_key = (file_base_name, scan_number)

            # Load data if not already loaded
            if current_file_data is None:
                current_file_name = file_list[current_file_index]
                with open(f"{self.mgf_path}/{current_file_name}") as f:
                    current_file_data = f.readlines()
                line_index = 0

            while True:
                # If the current file has reached the end
                if line_index >= len(current_file_data):
                    current_file_index += 1

                    # No more files to process
                    if current_file_index >= len(file_list):
                        logging.info("Reached the end of the files. No more data.")
                        break
                    current_file_name = file_list[current_file_index]
                    with open(f"{self.mgf_path}/{current_file_name}") as f:
                        current_file_data = f.readlines()
                    line_index = 0

                # Read the current line
                current_line = current_file_data[line_index].strip()

                if current_line.startswith('TITLE='):
                    # Extract scan number and file name from the TITLE
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
                        with open(f"{self.mgf_path}/{current_file_name}") as f:
                            current_file_data = f.readlines()
                        line_index = 0
                        continue

                    # Initialize the spectrum if the current spectrum matches
                    if spec_key == current_spec_key:
                        spectrum = self._initialize_spectrum(self.max_length)

                if spec_key == current_spec_key:
                    if current_line.startswith('END IONS'):
                        # Normalize and encode the spectrum and metadata
                        spectrum = self._normalize_spectrum(spectrum)
                        sequence_x = self._encode_sequence(peptide_x, charge)
                        sequence_y = self._encode_sequence(peptide_y, charge)
                        psm_features_x = self._add_psm_features(cscore_x, ifi_x, diff_rt_x, xcorr_x)
                        psm_features_y = self._add_psm_features(cscore_y, ifi_y, diff_rt_y, xcorr_y)
                        delta_features = self._add_delta_features(delta_cscore, delta_xcorr)  # count

                        # Yield the processed input data
                        yield {
                            input_spectrum: spectrum,
                            input_sequence_x: sequence_x,
                            input_sequence_y: sequence_y,
                            input_psm_features_x: psm_features_x,
                            input_psm_features_y: psm_features_y,
                            input_delta_features: delta_features
                        }

                        line_index += 1
                        break

                    # Collect peak data if it's a valid spectrum line
                    if current_line.strip() and current_line.split()[0].replace('.', '', 1).isdigit():
                        mz = float(current_line.split()[0])
                        intensity = float(current_line.split()[1])
                        loc = int(mz // self.mz_resolution)
                        if loc <= self.max_length:
                            spectrum[loc] += intensity

                # Move to the next line
                line_index += 1

    @staticmethod
    def _initialize_spectrum(length):
        # Initialize a zero-filled spectrum array of the given length
        return np.zeros(length, dtype=np.float32)

    @staticmethod
    def _normalize_spectrum(spectrum):
        # Normalize the spectrum by dividing by the maximum intensity
        max_intensity = np.max(spectrum)
        return spectrum / max_intensity if max_intensity > 0 else spectrum

    def _encode_sequence(self, sequence, charge):
        charge_idx = self.charge_index[int(charge)]

        # Pre-allocate a 40x28 matrix for sequence encoding
        encoded_sequence = np.zeros((40, 28), dtype=np.int8)

        for i, residue in enumerate(sequence):
            encoded_residue = encoded_sequence[i]
            # Mark amino acid position
            encoded_residue[self.amino_acid_index[residue]] = 1
            # Mark charge state position
            encoded_residue[charge_idx] = 1
            if residue in self.modification_index:
                # Mark modification position
                encoded_residue[self.modification_index[residue]] = 1

        return encoded_sequence

    @staticmethod
    def _add_psm_features(cscore, ifi, diff_rt, xcorr):
        # Combine PSM features into a single numpy array
        return np.array([cscore, ifi, diff_rt, xcorr], dtype=np.float32)

    @staticmethod
    def _add_delta_features(delta_cscore, delta_xcorr):
        # Combine delta features into a single numpy array
        return np.array([delta_cscore, delta_xcorr], dtype=np.float32)

    @staticmethod
    def _labeling(label):
        # Convert the label into a numpy array
        return np.array([label], dtype=np.int8)
