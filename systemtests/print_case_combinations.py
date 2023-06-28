import yaml
from pathlib import Path
from metadata_parser.Tutorial import Tutorials
from metadata_parser.Component import Components

tutorials_path= Path(__file__).parent.parent
systemtests_path = Path(__file__).parent

available_components = Components.from_yaml(systemtests_path / "components.yaml")

tutorial_dir = '../'
available_tutorials = Tutorials.from_path(tutorials_path)

tutorials = {}
for tutorial in available_tutorials:
    cases_combinations = [ f"{combination}" for combination in tutorial.case_combinations]
    tutorials[tutorial.path] = cases_combinations

print (yaml.dump(tutorials))