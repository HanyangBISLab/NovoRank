import os
import logging
import warnings
import tensorflow as tf

from tensorflow import keras
from src.model.preprocess import *

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class NovoRankTrainer:

    """
    A class to handle the training of the NovoRank model.
    """

    def __init__(self, train_data, val_data, configs):
        self.train_data = train_data
        self.val_data = val_data
        self.configs = configs

        self.gen_data_train = None
        self.gen_data_val = None
        self.model = None

        self.input_generator = InputDataGenerator(configs['path']['mgf_path'])

    def preprocess(self):

        """
        Preprocesses the training and validation data by selecting relevant columns.
        """

        selected_columns = ['Source File', 'Scan number', 'Peptide_x', 'Peptide_y',
                            'Charge_x', 'Score_x', 'Score_y', 'Delta Score_y',
                            'Normalized Internal Fragment Ions_x', 'Normalized Internal Fragment Ions_y',
                            'Difference_RT (min)_x', 'Difference_RT (min)_y',
                            'XCorr_x', 'XCorr_y', 'Delta XCorr', 'Label_x']

        # Filter the data to include only the selected columns
        self.train_data = self.train_data[selected_columns].copy()
        self.val_data = self.val_data[selected_columns].copy()

    def generator(self, batch_size):

        """
        Creates data generators for training and validation datasets.
        Converts the input data into a format suitable for training the model.
        """

        output_type = ({
            'input_spectrum': tf.float32, 'input_sequence_x': tf.float32, 'input_sequence_y': tf.float32,
            'input_psm_features_x': tf.float32, 'input_psm_features_y': tf.float32,
            'input_delta_features': tf.float32
        }, tf.float32)
        input_shape = ({
            'input_spectrum': [50000], 'input_sequence_x': [40, 28], 'input_sequence_y': [40, 28],
            'input_psm_features_x': [4], 'input_psm_features_y': [4], 'input_delta_features': [2]
        }, [1])

        # Training data generator
        self.gen_data_train = tf.data.Dataset.from_generator(
            lambda: self.input_generator.generate_data(self.train_data),
            output_types=output_type,
            output_shapes=input_shape
        ).batch(batch_size).prefetch(tf.data.experimental.AUTOTUNE)

        # Validation data generator
        self.gen_data_val = tf.data.Dataset.from_generator(
            lambda: self.input_generator.generate_data(self.val_data),
            output_types=output_type,
            output_shapes=input_shape
        ).batch(batch_size).prefetch(tf.data.experimental.AUTOTUNE)

    def combine_model(self, model_fn, model_fn_):

        """
        Generates the final model architecture.
        """

        # Define input layers
        spectra = keras.layers.Input(shape=(50000, 1), name='input_spectrum')
        sequence_x = keras.layers.Input(shape=(40, 28), name='input_sequence_x')
        sequence_y = keras.layers.Input(shape=(40, 28), name='input_sequence_y')
        psm_features_x = keras.layers.Input(shape=(4,), name='input_psm_features_x')
        psm_features_y = keras.layers.Input(shape=(4,), name='input_psm_features_y')
        delta_features = keras.layers.Input(shape=(2,), name='input_delta_features')

        base_model = model_fn
        model_x = base_model([spectra, sequence_x, psm_features_x])
        model_y = base_model([spectra, sequence_y, psm_features_y])

        output = model_fn_(model_x, model_y, delta_features)

        # Define the final model
        self.model = tf.keras.Model(
            inputs=[spectra, sequence_x, sequence_y, psm_features_x, psm_features_y, delta_features],
            outputs=output
        )

        # Compile the model with Adam optimizer and binary crossentropy loss
        self.model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

    def train_model(self, epoch):

        """
        Trains the model using the training and validation data generators.
        Optionally saves model checkpoints during training.
        """

        if self.configs['params']['checkpoint']:
            checkpoint_path = os.path.join(self.configs['path']['model_save']['path'],
                                           self.configs['path']['model_save']['filename'].split('.')[0]+'_{epoch}.h5')

            callback = tf.keras.callbacks.ModelCheckpoint(filepath=checkpoint_path,
                                                          save_weights_only=False, verbose=0, save_freq=1)

            # Train the model with checkpoint callback
            history = self.model.fit(
                self.gen_data_train,
                validation_data=self.gen_data_val,
                epochs=epoch,
                callbacks=[callback],
                verbose=1
            )
        else:
            # Train the model without checkpointing
            history = self.model.fit(
                self.gen_data_train,
                validation_data=self.gen_data_val,
                epochs=epoch,
                verbose=1
            )

        return history

    def save_model(self):

        """
        Saves the trained model to the specified path.
        """

        logger.info(f"Starting to save the model to {self.configs['path']['model_save']['path']}...")

        # Save the model
        self.model.save(os.path.join(self.configs['path']['model_save']['path'],
                                     self.configs['path']['model_save']['filename']))
        logger.info(f"Model save process completed successfully. The model has been saved as {self.configs['path']['model_save']['filename']}.")
