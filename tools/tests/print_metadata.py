from pathlib import Path
from metadata_parser.Tutorial import Tutorials
from metadata_parser.Component import Components
from systemtests.TestSuite import TestSuites

from paths import PRECICE_TUTORIAL_DIR

available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)
print("Fount the following tutorials read from the metadata.yaml")
for tutorial in available_tutorials:
    print(tutorial)