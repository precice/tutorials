import numpy as np
from typing import Dict


def switching_function(model_res: int, location: np.ndarray, t: float, input: Dict, prev_output: Dict):
    """
    Determine which model to use based on the resolution, location, time, input, and previous output.

    Parameters
    ----------
    model_res : int
        The resolution level of the model (0 for full-order model, >0 for reduced-order model).
    location : np.ndarray
        The macro-scale coordinates corresponding to the micro-scale simulation.
    t : float
        The current time step.
    input : Dict
        The input data for the model, including concentration.
    prev_output : dict
        The output from the previous time step.

    Returns
    -------
    model_res_offset : int
        The offset to change the model used from the provided hierarchy of models.
    """
    model_res_offset = 0

    # Only use the full-order model for the first time step
    if t == 0.0:
        if model_res > 0:
            model_res_offset = -1
        else:
            model_res_offset = 0

    # After the first time step, the model is selected dynamically based on
    # the concentration value and the model resolution.
    else:
        concentration = input["Concentration"]
        is_valid_range = 0.2 < concentration < 0.3
        is_fom = model_res == 0

        if is_fom and is_valid_range:
            model_res_offset = 1
        if not is_fom and not is_valid_range:
            model_res_offset = -1

    return model_res_offset
