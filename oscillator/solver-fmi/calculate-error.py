import numpy as np
import pandas as pd
import json
import os
import argparse
from enum import Enum
import sys


class Participant(Enum):
    MASS_LEFT = "Mass-Left"
    MASS_RIGHT = "Mass-Right"


parser = argparse.ArgumentParser()
parser.add_argument("fmi_setting_file_left", help="Path to the fmi setting file for MassLeft.", type=str)
parser.add_argument("precice_setting_file_left", help="Path to the precice setting file for MassLeft.", type=str)
parser.add_argument("fmi_setting_file_right", help="Path to the fmi setting file for MassRight.", type=str)
parser.add_argument("precice_setting_file_right", help="Path to the precice setting file for MassRight.", type=str)
parser.add_argument("participant_name", help="Participant for which the error should be calculated", type=str,
                    choices=[p.value for p in Participant])
args = parser.parse_args()

# Get input files
fmi_file_left = args.fmi_setting_file_left
precice_file_left = args.precice_setting_file_left
fmi_file_right = args.fmi_setting_file_right
precice_file_right = args.precice_setting_file_right
participant_name = args.participant_name

# Read json files
folder = os.path.dirname(os.path.join(os.getcwd(), os.path.dirname(sys.argv[0]), fmi_file_left))
path = os.path.join(folder, os.path.basename(fmi_file_left))
read_file = open(path, "r")
fmi_data_left = json.load(read_file)

folder = os.path.dirname(os.path.join(os.getcwd(), os.path.dirname(sys.argv[0]), precice_file_left))
path = os.path.join(folder, os.path.basename(precice_file_left))
read_file = open(path, "r")
precice_data_left = json.load(read_file)

folder = os.path.dirname(os.path.join(os.getcwd(), os.path.dirname(sys.argv[0]), fmi_file_right))
path = os.path.join(folder, os.path.basename(fmi_file_right))
read_file = open(path, "r")
fmi_data_right = json.load(read_file)

folder = os.path.dirname(os.path.join(os.getcwd(), os.path.dirname(sys.argv[0]), precice_file_right))
path = os.path.join(folder, os.path.basename(precice_file_right))
read_file = open(path, "r")
precice_data_right = json.load(read_file)


# Define variables
k_1 = fmi_data_left["model_params"]["spring_fixed.c"]
k_2 = fmi_data_right["model_params"]["spring_fixed.c"]
u0_1 = fmi_data_left["initial_conditions"]["mass.u"]
u0_2 = fmi_data_right["initial_conditions"]["mass.u"]
k_12_left = fmi_data_left["model_params"]["spring_middle.c"]
k_12_right = fmi_data_right["model_params"]["spring_middle.c"]

if k_12_left == k_12_right:
    k_12 = k_12_left
else:
    raise Exception("k_12 has to be equal in both participants. Please adjust input values.")


# Define analytical solution and read computed results
K = np.array([[k_1 + k_12, -k_12], [-k_12, k_2 + k_12]])
eigenvalues, eigenvectors = np.linalg.eig(K)
omega = np.sqrt(eigenvalues)
A, B = eigenvectors
c = np.linalg.solve(eigenvectors, [u0_1, u0_2])

if participant_name == Participant.MASS_LEFT.value:
    filename = fmi_data_left["simulation_params"]["output_file_name"]
    df = pd.read_csv(filename, delimiter=',')

    def u_analytical(t): return c[0] * A[0] * np.cos(omega[0] * t) + c[1] * A[1] * np.cos(omega[1] * t)

elif participant_name == Participant.MASS_RIGHT.value:
    filename = fmi_data_right["simulation_params"]["output_file_name"]
    df = pd.read_csv(filename, delimiter=',')

    def u_analytical(t): return c[0] * B[0] * np.cos(omega[0] * t) + c[1] * B[1] * np.cos(omega[1] * t)

times = df.iloc[:, 0]
positions = df.iloc[:, 1]


# Calculate error
error = np.max(abs(u_analytical(np.array(times)) - np.array(positions)))
print("Error w.r.t analytical solution:")
print(f"{error}")
