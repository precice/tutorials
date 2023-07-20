import yaml
from pathlib import Path
from metadata_parser.Tutorial import Tutorials
from metadata_parser.Component import Components
from paths import PRECICE_TUTORIAL_DIR,PRECICE_TESTS_DIR

available_components = Components.from_yaml(PRECICE_TESTS_DIR / "components.yaml")

available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)

tutorials = {}
for tutorial in available_tutorials:
    cases_combinations = [ f"{combination}" for combination in tutorial.case_combinations]
    tutorials[tutorial.path] = cases_combinations

print (yaml.dump(tutorials))