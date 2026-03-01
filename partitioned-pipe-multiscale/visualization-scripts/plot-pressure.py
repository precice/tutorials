#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# User options
# -----------------------------
xmin = None
xmax = None

# 3D centerline (OpenFOAM sets raw): 2 columns -> x p
FOAM_X_COL = 0
FOAM_P_COL = 1

# 1D final_fields.txt: columns -> x u p
TXT_X_COL = 0
TXT_P_COL = 2


FILES = {

    "1D-3D": (
        ("../fluid1d-left/final_fields.txt", "txt1d"),
        ("../fluid3d-right/postProcessing/sampleDict/5/centerline_p.xy", "foamXY"),
        "1D Left Pipe (0-20 m)",
        "3D Right Pipe (20-40 m)",
        "../images/pressure_distribution_1D3D.png",
    ),

    "3D-1D": (
        ("../fluid3d-left/postProcessing/sampleDict/5/centerline_p.xy", "foamXY"),
        ("../fluid1d-right/final_fields.txt", "txt1d"),
        "3D Left Pipe (0-20 m)",
        "1D Right Pipe (20-40 m)",
        "../images/pressure_distribution_3D1D.png",
    ),

    "1D-1D": (
        ("../fluid1d-left/final_fields.txt", "txt1d"),
        ("../fluid1d-right/final_fields.txt", "txt1d"),
        "1D Left Pipe (0-20 m)",
        "1D Right Pipe (20-40 m)",
        "../images/pressure_distribution_1D1D.png",
    ),

    "3D-3D": (
        ("../fluid3d-left/postProcessing/sampleDict/5/centerline_p.xy", "foamXY"),
        ("../fluid3d-right/postProcessing/sampleDict/5/centerline_p.xy", "foamXY"),
        "3D Left Pipe (0-20 m)",
        "3D Right Pipe (20-40 m)",
        "../images/pressure_distribution_3D3D.png",
    ),
}


PLOT_CASES = {
    "1D-3D",
    "3D-1D",
    # "1D-1D",
    # "3D-3D"
}


def crop(x, y):
    if xmin is None and xmax is None:
        return x, y
    m = np.ones_like(x, dtype=bool)
    if xmin is not None:
        m &= x >= xmin
    if xmax is not None:
        m &= x <= xmax
    return x[m], y[m]


def load_foam_xy(path):
    data = np.loadtxt(path, comments="#")
    x = data[:, FOAM_X_COL].astype(float)
    p = data[:, FOAM_P_COL].astype(float)
    o = np.argsort(x)
    return crop(x[o], p[o])


def load_txt1d(path):
    data = np.loadtxt(path, comments="#")
    x = data[:, TXT_X_COL].astype(float)
    p = data[:, TXT_P_COL].astype(float)
    o = np.argsort(x)
    return crop(x[o], p[o])


for case in PLOT_CASES:
    (pathL, kindL), (pathR, kindR), labelL, labelR, out = FILES[case]

    if kindL == "foamXY":
        xL, pL = load_foam_xy(pathL)
    else:
        xL, pL = load_txt1d(pathL)

    if kindR == "foamXY":
        xR, pR = load_foam_xy(pathR)
    else:
        xR, pR = load_txt1d(pathR)

    plt.figure(figsize=(8, 5))
    plt.plot(xL, pL, 'b-', label=labelL)
    plt.plot(xR, pR, 'r-', label=labelR)

    plt.xlabel('Axial coordinate [m]')
    plt.ylabel('Pressure [Pa]')
    plt.title(f'Pressure Distribution Along Coupled {case} Pipe')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()
