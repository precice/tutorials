"""
Micro simulation Surrogate, requrie previous computation of surrogate model
"""
import joblib
import numpy as np
import math


class MicroSimulation:

    def __init__(self, sim_id):
        """
        Get the micro-scale model from the BayesValidRox surrogate model.

        Parameters
        ----------
        sim_id : int
            The simulation ID for the micro-scale simulation.
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
        output_data["K00"] = 0.4912490635619572
        output_data["K11"] = 0.4912490635989945
        output_data["Porosity"] = 0.4933482661391027

        if self._sim_id == 0:
            output_data["K00"] = 0.4912490640081466
            output_data["K11"] = 0.4912490640081367

        return output_data

    def get_state(self):
        """
        Get the current state of the micro-scale simulation.

        Returns
        -------
        state : dict
            The current state of the micro-scale simulation.
        """
        return self._state

    def set_state(self, state):
        """
        Set the current state of the micro-scale simulation.

        Parameters
        ----------
        state : dict
            The state to set for the micro-scale simulation.
        """
        self._state = state

    def solve(self, macro_data, dt):
        """
        Solve the micro-scale simulation using the surrogate model.

        Parameters
        ----------
        macro_data : dict
            The macro-scale data required for the micro-scale simulation.
        dt : float
            The time step for the micro-scale simulation.

        Returns
        -------
        output_data : dict
            The output data from the micro-scale simulation.
        """
        model_eval, _ = self._model.eval_metamodel(np.array([macro_data["Concentration"]])[:, np.newaxis])
        output_data = dict()
        output_data["K00"] = model_eval["K00"][0][0]
        output_data["K11"] = model_eval["K11"][0][0]
        output_data["Porosity"] = model_eval["Porosity"][0][0]
        output_data["grain_size"] = math.sqrt((1 - model_eval["Porosity"][0][0]) / math.pi)

        return output_data
