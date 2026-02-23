from __future__ import division, print_function


def performAction(time, sourceData, targetData):
    timeThreshold = 0.2
    ramp = min(time / timeThreshold, 1.0)
    targetData[:] = ramp * sourceData
