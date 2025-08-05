---
title: Partitioned Water Hammer
keywords: OpenFOAM, python, preCICE, multiscale, fluid, FSI, transient
summary: The Partitioned Water Hammer tutorial simulates unsteady pressure wave propagation in pipe systems using different 1D and 3D configurations coupled via preCICE.
---

## Setup

This tutorial demonstrates how to model the **water hammer phenomenon** — a transient pressure surge caused by a rapid change in flow — using **partitioned coupling** with [preCICE](https://precice.org).

It supports multiple configurations:
- **1D–1D**: Python solvers on both sides
- **1D–3D**: Python 1D solver coupled to OpenFOAM
- **3D–3D**: OpenFOAM solvers on both sides

We exchange:
- **Velocity** from downstream to upstream
- **Pressure** from upstream to downstream

This allows realistic simulation of pressure wave propagation.

## Folder structure

```bash
water-hammer/
├── case-1d/            # Monolithic 1D simulation
|   └── fluid1d-python-uncoupled/
├── case-3d/            # Monolithic 3D simulation
│   └── fluid3d-openfoam-uncoupled/
├── case-1d-1d/         # Partitioned 1D–1D simulation
│   ├── fluid1dleft-python/
│   └── fluid1dright-python/
├── case-1d-3d/         # Partitioned 1D–3D simulation
│   ├── fluid1d-python/
│   └── fluid3d-openfoam/
├── case-3d-3d/         # Partitioned 3D–3D simulation
│   ├── fluid3d-openfoam-left/
│   └── fluid3d-openfoam-right/
├── results/       # Output data and plots
└── README.md      # This file
```

## How to use

### 1D–1D Coupling (Python–Python)

In two different terminals execute

```bash
cd case-1d-1d/fluid1dleft-python && ./run.sh
```

```bash
cd case-1d-1d/fluid1dright-python && ./run.sh
```

### 3D-3D Coupling (OpenFOAM-OpenFOAM)

In two different terminals execute

```bash
cd case-3d-3d/fluid3d-openfoam-left && ./run.sh
```

```bash
cd case-3d-3d/fluid3d-openfoam-right && ./run.sh
```

### 1D-3D Coupling (Python-OpenFOAM)

In two different terminals execute

```bash
cd case-1d-3d/fluid1d-python && ./run.sh
```

```bash
cd case-1d-3d/fluid3d-openfoam && ./run.sh
```

### 1D (Monolithic, Python)

In one terminal, execute

```bash
cd case-1d/fluid1d-python-uncoupled && ./run.sh
```

### 3D (Monolithic, OpenFOAM)

In one terminal, execute

```bash
cd case-3d/fluid3d-openfoam-uncoupled && ./run.sh
```