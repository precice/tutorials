#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# User options
# -----------------------------
tmin = None
tmax = None

FOAM_P_COL = 1
TXT_SEP = ";"
TXT_T_COL_1B = 1
TXT_P_COL_1B = 4

FILES = {
    "1D-3D": ("../fluid3d-right/postProcessing/probes/0/p", "foam"),
    "3D-1D": ("../fluid1d-right/probes.txt", "txt"),
    "1D-1D": ("../fluid1d-right/probes.txt", "txt"),
    "3D-3D": ("../fluid3d-right/postProcessing/probes/0/p", "foam"),
}

PLOT_CASES = {
    "1D-3D",
    "3D-1D",
    # "1D-1D",
    # "3D-1D"
}


def crop(t, y):
    if tmin is None and tmax is None:
        return t, y
    m = np.ones_like(t, dtype=bool)
    if tmin is not None:
        m &= t >= tmin
    if tmax is not None:
        m &= t <= tmax
    return t[m], y[m]


def dedup_sort_keep_last(t, y):
    # keep last value for duplicate times; then sort by time
    idx_last = {}
    for i, ti in enumerate(t):
        idx_last[ti] = i
    idx = np.fromiter(idx_last.values(), dtype=int)
    t2, y2 = t[idx], y[idx]
    o = np.argsort(t2)
    return t2[o], y2[o]


def load_foam(path, p_col=FOAM_P_COL):
    t, p = [], []
    for line in open(path, "r"):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        cols = s.split()
        t.append(float(cols[0]))
        p.append(float(cols[p_col]))
    t = np.asarray(t, float)
    p = np.asarray(p, float)
    t, p = dedup_sort_keep_last(t, p)
    return crop(t, p)


def load_txt(path, sep=TXT_SEP, t_col_1b=TXT_T_COL_1B, p_col_1b=TXT_P_COL_1B):
    ti, pi = t_col_1b - 1, p_col_1b - 1
    t, p = [], []
    for line in open(path, "r"):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        cols = [c.strip() for c in s.split(sep)]
        t.append(float(cols[ti]))
        p.append(float(cols[pi]))
    t = np.asarray(t, float)
    p = np.asarray(p, float)
    t, p = dedup_sort_keep_last(t, p)
    return crop(t, p)


for label in PLOT_CASES:
    path, kind = FILES[label]

    if kind == "foam":
        t, p = load_foam(path)
    else:
        t, p = load_txt(path)

    outlet_domain = "1D domain" if label.split("-")[1] == "1D" else "3D domain"
    out = f"../images/p_outlet_{label.replace('-', '')}.png"

    plt.figure(figsize=(9, 4.8))
    plt.plot(t, p, linewidth=1.5)
    plt.xlabel("Time [s]")
    plt.ylabel("Pressure p [Pa]")
    plt.title(f"Pressure evolution at the outlet of the {outlet_domain} ({label})")
    plt.grid(True, linewidth=0.3)
    plt.tight_layout()

    plt.savefig(out, dpi=200)
    plt.close()
