from pathlib import Path
from metadata_parser.Tutorial import Tutorials
from metadata_parser.Component import Components
from systemtests.TestSuite import TestSuites
tutorials_path= Path(__file__).parent.parent

available_tutorials = Tutorials.from_path(tutorials_path)
print("Fount the following tutorials read from the metadata.yaml")
for tutorial in available_tutorials:
    print(tutorial)