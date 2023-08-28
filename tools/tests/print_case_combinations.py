import yaml
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

    tutorials = {}
    for tutorial in available_tutorials:
        cases_combinations = [
            f"{combination}" for combination in tutorial.case_combinations]
        tutorials[tutorial.path.name] = cases_combinations

    print(yaml.dump(tutorials))


if __name__ == '__main__':
    main()
