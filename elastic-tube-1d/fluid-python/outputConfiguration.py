from __future__ import division
from enum import Enum


class PlottingModes(Enum):
    OFF = 0  # no plotting over time
    VIDEO = 1  # create a video
    DEBUG = 2  # provide a debug plot over time of the simulation


class OutputModes(Enum):
    OFF = 0  # no plotting over time
    VTK = 1  # produce VTK output
