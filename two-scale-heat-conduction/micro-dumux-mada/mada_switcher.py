import numpy as np


def switching_function(resolutions, locations, t, inputs, prev_output, active):
    output = np.zeros_like(resolutions)
    # in the beginning we only want FOM
    if t == 0.0:
        output[(resolutions > 0) * active] = -1
    # after small init phase, we want dynamic model selection
    # for test purposes we say ROM is accurate in range
    else:
        concentrations = np.array([d['concentration'] for d in inputs])
        valid_range_mask = (concentrations > 0.45) * (concentrations < 0.55)
        fom_mask = resolutions == 0

        output[fom_mask * valid_range_mask * active] = 1
        output[(~fom_mask) * (~valid_range_mask) * active] = -1

    return output
