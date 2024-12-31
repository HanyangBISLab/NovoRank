import tensorflow as tf

from tensorflow import keras
from tensorflow.keras.layers import Conv1D, LSTM, Bidirectional, Dense, Flatten  # BatchNormalization
from tensorflow.keras.layers import MaxPool1D, LeakyReLU, Dropout, Reshape, concatenate, Activation


class NovoRankModel:
    def __init__(self):

        """
        Defines the input shapes for the spectrum, sequence, and additional features.
        """

        self.spec_shape = (50000, 1)
        self.seq_shape = (40, 28)
        self.plus_shape = (4,)

    def build_model(self):

        """
        Builds the NovoRank model.
        Combines Conv1D and BiLSTM-based encoders with dense layers.
        """

        # Define input tensors
        encoded_spectrum = keras.layers.Input(shape=self.spec_shape)
        encoded_sequence = keras.layers.Input(shape=self.seq_shape)
        psm_features = keras.layers.Input(shape=self.plus_shape)

        # Encode spectrum and sequence using respective encoders
        conv_features = self.conv_encoder(encoded_spectrum)
        bilstm_features = self.bilstm_encoder(encoded_sequence)

        # Concatenate encoded features with additional PSM features
        output = concatenate([conv_features, bilstm_features, psm_features])

        # Add dense layers for feature fusion
        output = self.dense_block(output, 32)
        output = self.dense_block(output, 32)
        output = Dense(16)(output)

        # Reshape the output to flatten it
        output = Reshape((-1,))(output)

        # Create the model
        model = tf.keras.Model(inputs=[encoded_spectrum, encoded_sequence, psm_features], outputs=output)

        return model

    def conv_encoder(self, spectrum_input):

        """
        Builds a Conv1D-based encoder for spectrum data.
        """

        # First convolutional block
        conv = Conv1D(8, 30, strides=1, padding="valid")(spectrum_input)
        conv = LeakyReLU(alpha=0.01)(conv)
        conv = Dropout(0.1)(conv)
        conv = MaxPool1D(pool_size=30, strides=30, padding='valid')(conv)

        # Second convolutional block
        conv = Conv1D(16, 30, strides=1, padding="valid")(conv)
        conv = LeakyReLU(alpha=0.01)(conv)
        conv = Dropout(0.1)(conv)
        conv = MaxPool1D(pool_size=30, strides=30, padding='valid')(conv)

        # Flatten and pass through dense blocks
        conv = Flatten()(conv)
        conv = self.dense_block(conv, 16)
        conv = self.dense_block(conv, 16)
        conv = Dense(16)(conv)

        return conv

    def bilstm_encoder(self, sequence_input):

        """
        Builds a Bidirectional LSTM-based encoder for sequence data.
        """

        # BiLSTM layer
        bilstm = Bidirectional(LSTM(8), merge_mode='concat')(sequence_input)

        # Pass through dense blocks
        bilstm = self.dense_block(bilstm, 16)
        bilstm = self.dense_block(bilstm, 16)
        bilstm = Dense(16)(bilstm)

        return bilstm

    @staticmethod
    def dense_block(input_layer, units, dropout_rate=0.1, leaky_relu_alpha=0.01):

        """
        Creates a dense block with a LeakyReLU activation and dropout.
        """

        x = Dense(units)(input_layer)
        x = LeakyReLU(alpha=leaky_relu_alpha)(x)
        x = Dropout(dropout_rate)(x)

        return x

    def build_output_model(self, model_x, model_y, delta_features):

        """
        Builds output model by combining two models' outputs with delta features.
        """

        # Concatenate all features
        features = concatenate([model_x, model_y, delta_features])

        # Add dense layers for final processing
        output = self.dense_block(features, 32)
        output = self.dense_block(output, 32)
        output = Dense(1)(output)
        output = Activation('sigmoid')(output)

        return output
