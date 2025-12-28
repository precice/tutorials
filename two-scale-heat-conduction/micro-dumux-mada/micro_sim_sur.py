"""
Micro simulation Surrogate, requrie previous computation of surrogate model
"""
import os
import subprocess
from bayesvalidrox import PyLinkForwardModel, Input, PCE, ExpDesigns, Engine
import h5py
import joblib
import numpy as np
import math


class MicroSimulation:

    def __init__(self, sim_id):
        """
        Constructor of MicroSimulation class.
        """
        self._sim_id = sim_id
        self._state = None

        self._model = None
        with open('micro-dumux-surrogate.pkl', 'rb') as input:
            self._model = joblib.load(input)
        if self._model is None:
            raise RuntimeError("Failed to load model.")

    def initialize(self):
        output_data = dict()
        output_data["k_00"] = 0.4912490635619572
        output_data["k_11"] = 0.4912490635989945
        output_data["porosity"] = 0.4933482661391027

        if self._sim_id == 0:
            output_data["k_00"] = 0.4912490640081466
            output_data["k_11"] = 0.4912490640081367

        return output_data

    def get_state(self):
        return self._state

    def set_state(self, state):
        self._state = state

    def solve(self, macro_data, dt):
        model_eval, _ = self._model.eval_metamodel(np.array([macro_data["concentration"]])[:, np.newaxis])
        output_data = dict()
        output_data["k_00"] = model_eval["k_00"][0][0]
        output_data["k_01"] = model_eval["k_01"][0][0]
        output_data["k_10"] = model_eval["k_10"][0][0]
        output_data["k_11"] = model_eval["k_11"][0][0]
        output_data["porosity"] = model_eval["porosity"][0][0]
        output_data["grain_size"] = math.sqrt((1 - model_eval["porosity"][0][0]) / math.pi)

        return output_data
