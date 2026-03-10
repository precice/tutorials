#!/usr/bin/env python
import sys
from pathlib import Path
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


def resolve_input(case_dir):
    parts = case_dir.split("-")
    last = parts[-1]

    if last == "openfoam":
        data_file = Path(case_dir) / "postProcessing" / "probes" / "0" / "p"
        loader = load_foam
    elif last == "nutils":
        data_file = Path(case_dir) / "probes.txt"
        loader = load_txt
    else:
        raise ValueError(f"Unknown case type '{last}'. Expected a directory ending in ...-openfoam' or '...-nutils'.")

    return data_file, loader


if len(sys.argv) == 1:
    case_dir = "fluid3d-right-openfoam"
    case_name = case_dir
    print(f"[INFO] No directory provided. Using default: {case_dir}")
elif len(sys.argv) == 2:
    case_dir = Path(sys.argv[1])
    case_name = case_dir.name
else:
    print("Usage: python plot-pressure.py <case-directory>")
    print("Example: python plot-pressure.py fluid3d-right-openfoam")
    sys.exit(1)

allowed_cases = ["fluid3d-right-openfoam", "fluid1d-right-nutils"]

if case_name not in allowed_cases:
    print(f"[ERROR] Invalid case directory: {case_dir}")
    print("Allowed cases:")
    for c in allowed_cases:
        print(f"  - {c}")
    sys.exit(1)

try:
    data_file, loader = resolve_input(case_name)
except ValueError as e:
    print(f"[ERROR] {e}")
    sys.exit(1)

if not data_file.is_file():
    print(f"[ERROR] File not found: {data_file}")
    sys.exit(1)

t, p = loader(data_file)
out = f"images/p_outlet_{case_name}.png"

plt.figure(figsize=(9, 4.8))
plt.plot(t, p, linewidth=1.5)
plt.xlabel("Time [s]")
plt.ylabel("Pressure p [Pa]")
plt.title(f"Pressure evolution at the outlet ({case_name})")
plt.grid(True, linewidth=0.3)
plt.tight_layout()

plt.savefig(out, dpi=200)
plt.close()
print(f"[INFO] Image saved to: {out}")
