from metadata_parser.metdata import Tutorials
from paths import PRECICE_TUTORIAL_DIR

import argparse
import logging


def main():
    parser = argparse.ArgumentParser(description='Prints available Metadata for tutorials')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO', help='Set the logging level')
    args = parser.parse_args()

    # Configure logging based on the provided log level
    logging.basicConfig(level=args.log_level, format='%(levelname)s: %(message)s')

    print(f"Using log-level: {args.log_level}")

    available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)
    print("Fount the following tutorials read from the metadata.yaml")
    for tutorial in available_tutorials:
        print(tutorial)


if __name__ == '__main__':
    main()
