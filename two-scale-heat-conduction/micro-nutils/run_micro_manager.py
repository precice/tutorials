"""
Script to run the Micro Manager
"""

from micro_manager import MicroManager

manager = MicroManager("micro-manager-config.json")

manager.solve()
