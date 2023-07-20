import yaml
from pathlib import Path
from metadata_parser.Tutorial import Tutorials
from metadata_parser.Component import Components
from systemtests.TestSuite import TestSuites
tutorials_path= Path(__file__).parent.parent
systemtests_path = Path(__file__).parent

available_tutorials = Tutorials.from_path(tutorials_path)
available_components = Components.from_yaml(systemtests_path / "components.yaml")
available_testsuites = TestSuites.from_yaml(systemtests_path / "tests.yaml",available_tutorials)





print (available_testsuites)