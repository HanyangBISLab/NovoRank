import os
import warnings

from src.model.train import *
from src.model.inference import *
from src.model.preprocess import *
from src.model.base_model import *
from src.loader.dataloader import *
from src.utils import config_second, utils

# Ignore: warning, info massage
# Print: error massage
# os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

# Ignore: ALL massage
# Suppress TensorFlow warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Hide DeprecationWarnings occurring in Keras
warnings.filterwarnings("ignore", category=DeprecationWarning)

# warnings.filterwarnings('ignore')

# Future Tasks:
# Specify input and output types for the function
if __name__ == '__main__':
    args = config_second.parse_arguments()
    config_path = args.config
    configs = config_second.load_config(config_path)
    logging.info("Configuration loaded successfully.")

    # # Define classes
    interim_report_file_handler = FileHandler(os.path.join(configs['path']['save_path'],
                                                           configs['path']['interim_report']))
    de_novo_file_handler = FileHandler(os.path.join(configs['path']['search_results']['de_novo']['path'],
                                                    configs['path']['search_results']['de_novo']['filename']))
    processor = CrossCorrelationResultProcessor(configs['path']['xcorr_results_path'])
    preparation = DataPreparation()

    logging.info(f"TensorFlow version: {tf.__version__}")
    logging.info(f"Is TensorFlow built with CUDA: {tf.test.is_built_with_cuda()}")
    logging.info(f"TensorFlow build info: {tf.sysconfig.get_build_info()}")

    gpus = tf.config.experimental.list_physical_devices('GPU')
    logging.info(f"Available GPUs: {gpus}")

    if gpus:
        try:
            os.environ["CUDA_VISIBLE_DEVICES"] = "0"
            tf.config.experimental.set_memory_growth(gpus[0], True)
            logging.info("GPU configuration completed.")
        except RuntimeError as e:
            logging.error(f"Error configuring GPU: {e}")
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
        logging.info("No GPU found, using CPU.")

    logging.info("Starting dataset preparation process...")
    # Load the dataset and extract cross-correlation information
    dataset = interim_report_file_handler.load_csv()
    xcorr_df = processor.extract_cross_correlation_info()

    """
    Main function. Executes the training or inference process based on the TRAIN flag.
    """

    val_size = args.val_size
    batch_size = args.batch_size
    epoch = args.epoch

    # Process based on the TRAIN flag
    if configs['params']['train']:
        logging.info("Starting training process...")
        merged_df = preparation.get_train_dataset(dataset, xcorr_df)

        test_set = configs['params']['test_set']
        if test_set:
            logging.info("Creating train, validation, and test datasets...")
            train_df, val_df, test_df = preparation.train_val_split(merged_df, val_size, test_set)

            # Save datasets to CSV files
            train_df.to_csv(os.path.join(configs['path']['save_path'], 'train_data.csv'), index=False)
            val_df.to_csv(os.path.join(configs['path']['save_path'], 'val_data.csv'), index=False)
            test_df.to_csv(os.path.join(configs['path']['save_path'], 'test_data.csv'), index=False)
            logging.info("Train, validation, and test datasets created and saved successfully.")
        else:
            logging.info("Creating train and validation datasets...")
            train_df, val_df = preparation.train_val_split(merged_df, val_size)

            # Save datasets to CSV files
            train_df.to_csv(os.path.join(configs['path']['save_path'], 'train_data.csv'), index=False)
            val_df.to_csv(os.path.join(configs['path']['save_path'], 'val_data.csv'), index=False)
            logging.info("Train and validation datasets created and saved successfully.")

        logging.info("Dataset preparation process completed successfully.")

        # Model definition
        logging.info("Building model...")
        novorank = NovoRankModel()
        base_model = novorank.build_model()
        output_model = novorank.build_output_model
        logging.info("Model built successfully.")

        # Training process definition
        logging.info("Starting training preparation...")
        trainer = NovoRankTrainer(train_df, val_df, configs)
        trainer.preprocess()
        trainer.generator(batch_size)
        trainer.combine_model(base_model, output_model)
        logging.info("Training preparation completed successfully.")

        # Model training
        logging.info("Starting model training...")
        history = trainer.train_model(epoch)
        logging.info("Model training completed successfully.")

        # Model saving
        trainer.save_model()
    else:
        # Load the de novo dataset
        de_novo = de_novo_file_handler.load_csv()

        # Initiating the inference process
        logging.info("Starting inference process...")
        inference_df, above_max_seq_df, missing_xcorr_df = preparation.get_inference_dataset(dataset, xcorr_df)

        # Inference process definition
        logging.info("Starting inference preparation...")
        inference = NovoRankInference(inference_df, configs)
        inference.preprocess()
        inference.generator(batch_size)
        logging.info("Inference preparation completed successfully.")

        # Inference
        logging.info("Starting inference...")
        inference.run_inference(utils.top1(de_novo), above_max_seq_df, missing_xcorr_df)
        logging.info("Inference completed successfully.")

        # results saving
        # top 10 results?
        inference.save_results()
