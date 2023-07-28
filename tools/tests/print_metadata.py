from metadata_parser.metdata import Tutorials
from paths import PRECICE_TUTORIAL_DIR

available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)
print("Fount the following tutorials read from the metadata.yaml")
for tutorial in available_tutorials:
    print(tutorial)