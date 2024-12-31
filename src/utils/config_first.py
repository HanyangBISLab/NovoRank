import os
import yaml
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate candidates and extract features")
    parser.add_argument(
        '--search_ppm',
        '-search_ppm',
        type=float,
        required=True,
        help='Precursor tolerance (ppm) used in the de novo search.'
    )

    parser.add_argument(
        '--elution_time',
        '-elution_time',
        type=float,
        required=True,
        help='A total elution time (minutes) in the mass spectrometry assay.'
    )

    parser.add_argument(
        '--config',
        '-config',
        type=str,
        default='./config.yaml',
        help='Path to the config YAML file. (default: ./config.yaml)'
    )

    parser.add_argument(
        '--cluster_rt',
        '-cluster_rt',
        type=float,
        default=2,
        help='Retention time (in minutes) param for refining the clustering process. Default is 2 minutes.'
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
