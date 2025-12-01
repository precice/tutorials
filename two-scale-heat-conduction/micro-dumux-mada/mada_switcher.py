import numpy as np


def switching_function(resolution, location, t, input, prev_output):
    result = 0
    # in the beginning we only want FOM
    if t == 0.0:
        if resolution > 0: result = -1
        else: result = 0
    # after small init phase, we want dynamic model selection
    # for test purposes we say ROM is accurate in range
    else:
        concentration = input['concentration']
        is_valid_range = 0.45 < concentration < 0.55
        is_fom = resolution == 0

        if is_fom and is_valid_range: result = 1
        if not is_fom and not is_valid_range: result = -1

    return result
