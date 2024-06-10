"""
Script to run the Micro Manager
"""

from micro_manager import MicroManager
from argparse import ArgumentParser

print("Entered the Python scrip to run the Micro Manager")

parser = ArgumentParser()
parser.add_argument("--config", help="Path to the micro manager configuration file")
args = parser.parse_args()

print("Accepted the config file")

manager = MicroManager(args.config)

print("Micro Manager object created")

manager.initialize()

print("Micro Manager initialized")

manager.solve()
