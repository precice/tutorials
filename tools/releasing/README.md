# Releasing helpers

## Python dependency reference files

Each Python participant lists dependencies in its local `requirements.txt` using loose version ranges (for example `numpy >1, <2`, `pyprecice~=3.0`). Tutorial `run.sh` scripts install from these files as before.

Next to each `requirements.txt`, a sibling `requirements-reference.txt` records the latest PyPI versions that satisfy those constraints at generation time. Use it for reproducible installs:

```bash
pip install -r requirements-reference.txt
```

This tool lives under `tools/releasing/` and depends on:

```bash
pip install -r tools/releasing/requirements.txt
```

Update reference files after changing `requirements.txt` files. By default, only files whose resolved package pins changed are rewritten (timestamps alone do not create a diff). Use `--all` to force a full refresh, or pass a path to limit the update to one directory:

```bash
python3 tools/releasing/update-requirements-reference.py
python3 tools/releasing/update-requirements-reference.py path/to/participant
python3 tools/releasing/update-requirements-reference.py --all
python3 tools/releasing/update-requirements-reference.py --check
python3 tools/releasing/update-requirements-reference.py --check --fail-on-outdated
```

`--check` fails when a sibling `requirements-reference.txt` is missing, and warns when pins are outdated. `--fail-on-outdated` turns outdated pins into errors (used for release PRs to `master`).

Legacy FEniCS packages (for example `fenics-dolfin`) are installed from the system or PPA in several tutorials and appear in the reference files as non-PyPI dependencies.

When preparing a release pull request to `master`, regenerate and commit the sibling `requirements-reference.txt` files so the pinned snapshots match the current tutorial constraints.
