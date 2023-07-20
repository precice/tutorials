import yaml
from pathlib import Path
from metadata_parser.Tutorial import Tutorials
from metadata_parser.Component import Components
from systemtests.TestSuite import TestSuites
from paths import PRECICE_TESTS_DIR, PRECICE_TUTORIAL_DIR

available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)
available_components = Components.from_yaml(PRECICE_TESTS_DIR / "components.yaml")
available_testsuites = TestSuites.from_yaml(PRECICE_TESTS_DIR / "tests.yaml",available_tutorials)





print (available_testsuites)