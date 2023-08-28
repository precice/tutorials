from metadata_parser.metdata import Tutorials
from systemtests.TestSuite import TestSuites
from paths import PRECICE_TESTS_DIR, PRECICE_TUTORIAL_DIR


import argparse
import logging


def main():
    parser = argparse.ArgumentParser(description='Prints available Test Suites')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO', help='Set the logging level')
    args = parser.parse_args()

    # Configure logging based on the provided log level
    logging.basicConfig(level=args.log_level, format='%(levelname)s: %(message)s')

    print(f"Using log-level: {args.log_level}")

    available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)
    available_testsuites = TestSuites.from_yaml(
        PRECICE_TESTS_DIR / "tests.yaml", available_tutorials)

    print(available_testsuites)


if __name__ == '__main__':
    main()
