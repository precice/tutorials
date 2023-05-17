import glob
import os
import yaml
import itertools
import argparse


from classes import Tutorials , Components

available_components = Components.from_yaml("./components.yaml")

tutorial_dir = '../'
available_tutorials = Tutorials.from_path("../")

print("The following tutorials are described by the metadata")
for tutorial in available_tutorials.tutorials:
    print(tutorial)