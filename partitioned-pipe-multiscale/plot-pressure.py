#!/usr/bin/env python
import sys
import os
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
        ("fluid1d-left-nutils/final_fields.txt", "txt1d"),
        ("fluid3d-right-openfoam/postProcessing/sampleDict/5/centerline_p.xy", "foamXY"),
        "1D Left Pipe (0-20 m)",
        "3D Right Pipe (20-40 m)",
        "images/pressure_distribution_1D3D.png",
    ),

    "3D-1D": (
        ("fluid3d-left-openfoam/postProcessing/sampleDict/5/centerline_p.xy", "foamXY"),
        ("fluid1d-right-nutils/final_fields.txt", "txt1d"),
        "3D Left Pipe (0-20 m)",
        "1D Right Pipe (20-40 m)",
        "images/pressure_distribution_3D1D.png",
    ),

    "1D-1D": (
        ("fluid1d-left-nutils/final_fields.txt", "txt1d"),
        ("fluid1d-right-nutils/final_fields.txt", "txt1d"),
        "1D Left Pipe (0-20 m)",
        "1D Right Pipe (20-40 m)",
        "images/pressure_distribution_1D1D.png",
    ),

    "3D-3D": (
        ("fluid3d-left-openfoam/postProcessing/sampleDict/5/centerline_p.xy", "foamXY"),
        ("fluid3d-right-openfoam/postProcessing/sampleDict/5/centerline_p.xy", "foamXY"),
        "3D Left Pipe (0-20 m)",
        "3D Right Pipe (20-40 m)",
        "images/pressure_distribution_3D3D.png",
    ),
}


ARG_TO_CASE = {
    "1d3d": "1D-3D",
    "3d1d": "3D-1D",
    "1d1d": "1D-1D",
    "3d3d": "3D-3D",
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


if len(sys.argv) == 1:
    case_arg = "1d3d"
    print(f"[INFO] No case provided. Using default: {case_arg}")
elif len(sys.argv) == 2:
    case_arg = sys.argv[1].lower()
else:
    print("Usage: python plot-pressure.py <case>")
    print("Allowed cases: 1d3d, 3d1d, 1d1d, 3d3d")
    sys.exit(1)

if case_arg not in ARG_TO_CASE:
    print(f"[ERROR] Invalid case: {case_arg}")
    print("Allowed cases: 1d3d, 3d1d, 1d1d, 3d3d")
    sys.exit(1)

case = ARG_TO_CASE[case_arg]

(pathL, kindL), (pathR, kindR), labelL, labelR, out = FILES[case]

if not os.path.isfile(pathL):
    print(f"[ERROR] File not found: {pathL}")
    sys.exit(1)

if not os.path.isfile(pathR):
    print(f"[ERROR] File not found: {pathR}")
    sys.exit(1)

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

print(f"[INFO] Image saved to: {out}")
