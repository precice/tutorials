
import argparse
from pathlib import Path
from metadata_parser.metdata import Tutorials,Components
from systemtests.TestSuite import TestSuites




from paths import PRECICE_TUTORIAL_DIR,PRECICE_TESTS_RUN_DIR,PRECICE_TESTS_DIR

parser = argparse.ArgumentParser(description='generate reference data')

# Parse the command-line arguments
args = parser.parse_args()

available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)

available_components = Components.from_yaml(PRECICE_TESTS_DIR / "components.yaml")
test_suites = TestSuites.from_yaml(PRECICE_TESTS_DIR / "tests.yaml", available_tutorials)
print(test_suites)
