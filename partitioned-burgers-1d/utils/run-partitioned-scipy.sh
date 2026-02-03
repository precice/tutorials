#!/usr/bin/env bash
set -e -u

cd "$(dirname "$0")/.."

./clean-tutorial.sh

solver_path="./solver-scipy"

if [ -d "$solver_path/.venv" ]; then
	echo "Using existing virtual environment"
	. "$solver_path/.venv/bin/activate"
else
	echo "Creating new virtual environment"
	python3 -m venv "$solver_path/.venv"
	. "$solver_path/.venv/bin/activate"
	pip install -r $solver_path/requirements.txt && pip freeze > $solver_path/pip-installed-packages.log
fi

python3 utils/generate_ic.py --epoch "${1:-0}"

# full domain reference solution
echo "Running monolithic reference solution..."
(cd solver-scipy && ./run.sh)

(cd dirichlet-scipy && ./run.sh) &
(cd neumann-scipy && ./run.sh)

python3 utils/visualize_partitioned_domain.py