import os
import yaml
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate candidates and extract features")
    parser.add_argument(
        '--config',
        '-config',
        type=str,
        default='./config.yaml',
        help='Path to the config YAML file. (default: ./config.yaml)'
    )

    parser.add_argument(
        '--val_size',
        '-val_size',
        type=float,
        default=0.1,
        help='Validation set ratio used during training. Specify a value between 0 and 1. Default is 0.1.'
    )

    parser.add_argument(
        '--batch_size',
        '-batch_size',
        type=int,
        default=64,
        help='Batch size for dataset processing. Default is 64.'
    )

    parser.add_argument(
        '--epoch',
        '-epoch',
        type=int,
        default=10,
        help='Number of training epochs. Default is 10.'
    )

    return parser.parse_args()


def load_config(config_path):
    # Check if the path exists
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Path does not exist: {config_path}")

    # Check if the path is a file
    elif not os.path.isfile(config_path):
        raise TypeError(f"Path is not a file: {config_path}")

    # Load the configuration
    with open(config_path, "r", encoding='utf-8') as file:
        return yaml.safe_load(file)
