import os
import logging
import numpy as np
import pandas as pd
import tensorflow as tf

from src.model.preprocess import *

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class NovoRankInference:

    """
    A class to handle the inference of the NovoRank model.
    """

    def __init__(self, inference_data, configs):
        self.inference_data = inference_data
        self.configs = configs

        self.gen_data_inference = None
        # Load the pretrained model
        self.model = tf.keras.models.load_model(
            os.path.join(configs['path']['pretrained_model']['path'],
                         configs['path']['pretrained_model']['filename']),
            compile=False
        )

        self.input_generator = InputDataGenerator(configs['path']['mgf_path'])
        self.results_final = None

    def preprocess(self):

        """
        Preprocesses the inference data by selecting relevant columns.
        """

        selected_columns = ['Source File', 'Scan number', 'Peptide_x', 'Peptide_y',
                            'Charge_x', 'Score_x', 'Score_y', 'Delta Score_y',
                            'Normalized Internal Fragment Ions_x', 'Normalized Internal Fragment Ions_y',
                            'Difference_RT (min)_x', 'Difference_RT (min)_y',
                            'XCorr_x', 'XCorr_y', 'Delta XCorr']

        # Filter the data to include only the selected columns
        self.inference_data = self.inference_data[selected_columns].copy()

    def generator(self, batch_size):

        """
        Creates a data generator for the inference dataset.
        Converts the input data into a format suitable for inference.
        """

        input_names = [input_layer.name for input_layer in self.model.inputs]
        dimensions = [[50000], [40, 28], [40, 28], [4], [4], [2]]

        output_type = {name: tf.float32 for name in input_names}
        input_shape = {input_names[i]: dimensions[i] for i in range(len(input_names))}

        # Inference data generator
        self.gen_data_inference = tf.data.Dataset.from_generator(
            lambda: self.input_generator.generate_inference_data(self.inference_data, input_names),
            output_types=output_type,
            output_shapes=input_shape
        ).batch(batch_size).prefetch(tf.data.experimental.AUTOTUNE)

    def run_inference(self, de_novo, above_max_seq_df, missing_xcorr_df):

        """
        Runs the inference on the data using the pretrained model and processes the results.
        """

        # Set the model to inference mode (non-trainable)
        self.model.trainable = False
        # Predict using the inference data generator
        predictions = self.model.predict(self.gen_data_inference, verbose=1)

        # Add prediction scores to the data
        self.inference_data['Prediction Score'] = predictions
        self.inference_data['Prediction Binary'] = self._convert_predictions_to_binary(predictions)
        self.inference_data['Selected Peptide'] = self._select_peptide_by_prediction(self.inference_data)
        self.inference_data['Classification Score'] = self._calculate_cls_score(self.inference_data)

        # Combine results from different sources (missing XCorr, above max seq, and inference data)
        missing_xcorr_df = missing_xcorr_df[['Source File', 'Scan number', 'Peptide', 'Score']]
        above_max_seq_df = above_max_seq_df[['Source File', 'Scan number', 'Peptide_x', 'Score_x']].rename(
            columns={'Peptide_x': 'Peptide', 'Score_x': 'Score'})
        inference_data = self.inference_data[['Source File', 'Scan number',
                                              'Selected Peptide', 'Classification Score']].rename(
            columns={'Selected Peptide': 'Peptide', 'Classification Score': 'Score'})

        # Concatenate all results and sort by 'Source File' and 'Scan number'
        results_temp = pd.concat([missing_xcorr_df, above_max_seq_df, inference_data]).sort_values(
            by=['Source File', 'Scan number']).reset_index(drop=True)

        # Re-ranking the results based on the de novo data and threshold
        self.results_final = self._re_ranking(de_novo, results_temp, threshold=40)

    @staticmethod
    def _convert_predictions_to_binary(predictions):

        """
        Converts prediction scores into binary values (0 or 1) based on a threshold of 0.5.
        """

        return (predictions > 0.5).astype(int)

    @staticmethod
    def _select_peptide_by_prediction(dataset):

        """
        Selects the appropriate peptide based on the binary prediction.
        """

        return np.where(dataset['Prediction Binary'] == 0, dataset['Peptide_x'], dataset['Peptide_y'])

    @staticmethod
    def _calculate_cls_score(dataset):

        """
        Calculates the classification score based on the prediction score and binary prediction.
        """

        # Vectorized calculation using NumPy for efficient performance
        scores = np.where(
            dataset['Prediction Binary'] == 0,
            2 * (0.5 - dataset['Prediction Score']) * 100,
            2 * (dataset['Prediction Score'] - 0.5) * 100
        )

        return scores

    @staticmethod
    def _re_ranking(original, ranking, threshold=40):

        """
        Re-ranking the results based on a threshold
        If the ranking score exceeds the threshold, the peptide from the ranking is selected;
        otherwise, the peptide from the original data is kept.
        """

        # Merge the original data with the ranking data
        merged_data = pd.merge(original, ranking, on=['Source File', 'Scan number'], how='outer')

        # Extract relevant columns from the merged data and rename them for clarity
        ranking_data = merged_data[['Source File', 'Scan number', 'Peptide_x', 'Peptide_y',
                                    'Score_x', 'Score_y']].rename(
            columns={'Score_x': 'Score_o', 'Score_y': 'Score_r'})

        # Filter rows where both peptides are available
        ranking_data = ranking_data[ranking_data['Peptide_x'].notnull() & ranking_data['Peptide_y'].notnull()]

        # Apply the threshold: select peptide from ranking if the ranked score is above the threshold,
        # otherwise keep the original peptide
        ranking_data['Peptide'] = ranking_data.apply(
            lambda ser: ser['Peptide_y'] if ser['Score_r'] >= threshold else ser['Peptide_x'], axis=1)

        final_ranking = ranking_data[['Source File', 'Scan number', 'Peptide', 'Score_o']].rename(
            columns={'Score_o': 'Score'})

        # Handle cases where only one peptide is available (either from the original or the ranking)
        original_only = ranking_data[ranking_data['Peptide_x'].notnull() & ranking_data['Peptide_y'].isnull()]
        original_only = original_only[['Source File', 'Scan number', 'Peptide_x', 'Score_o']].rename(
            columns={'Peptide_x': 'Peptide', 'Score_o': 'Score'})

        ranking_only = ranking_data[ranking_data['Peptide_x'].isnull() & ranking_data['Peptide_y'].notnull()]
        ranking_only = ranking_only[['Source File', 'Scan number', 'Peptide_y', 'Score_r']].rename(
            columns={'Peptide_y': 'Peptide', 'Score_r': 'Score'})

        # Concatenate all the results and sort them by 'Source File' and 'Scan number'
        results = pd.concat([final_ranking, original_only, ranking_only]).sort_values(
            by=['Source File', 'Scan number']).reset_index(drop=True)

        return results

    def save_results(self):

        """
        Saves the inference results to a CSV file.
        """

        if self.results_final is None:
            raise ValueError("No results to save. Please run inference first.")

        # Define save path and file name
        save_path = self.configs['path']['save_path']
        save_name = self.configs['path']['final_report']

        # Save the final results to a CSV file
        logger.info(f"Starting to save results to {save_path}...")
        self.results_final.to_csv(os.path.join(save_path, save_name), index=False)
        logger.info(f"Results save process completed successfully. The results have been saved as {save_name}.")
