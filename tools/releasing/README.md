# Releasing helpers

## Python dependency reference files

Each Python participant lists dependencies in its local `requirements.txt` using loose version ranges (for example `numpy >1, <2`, `pyprecice~=3.0`). Tutorial `run.sh` scripts install from these files as before.

Next to each `requirements.txt`, a sibling `requirements-reference.txt` records the latest PyPI versions that satisfy those constraints at generation time. Use it for reproducible installs:

```bash
pip install -r requirements-reference.txt
```

Regenerate all reference files after changing tutorial requirements:

```bash
python3 tools/releasing/report_tutorial_requirements.py
python3 tools/releasing/report_tutorial_requirements.py --check   # CI-style check
```

Legacy FEniCS packages (for example `fenics-dolfin`) are installed from the system or PPA in several tutorials and appear in the reference files as non-PyPI dependencies.

When preparing a release pull request to `master`, regenerate and commit the sibling `requirements-reference.txt` files so the pinned snapshots match the current tutorial constraints.
