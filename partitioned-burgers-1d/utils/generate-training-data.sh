#!/usr/bin/env bash
set -e -u

# Execute this script from the tutorial root directory

mkdir -p solver-scipy/data-training

if [ ! -v PRECICE_TUTORIALS_NO_VENV ]
then
    if [ ! -d ".venv" ]; then
        python3 -m venv .venv
        source .venv/bin/activate
        pip install -r utils/requirements.txt && pip freeze > pip-installed-packages.log
    else
        source .venv/bin/activate
    fi
fi

# Number of training runs to generate
NUM_RUNS=200

echo "Generating ${NUM_RUNS} training data samples..."

for i in $(seq 0 $((NUM_RUNS - 1))); do
  echo "--- Generating epoch ${i} ---"

  # Generate IC
  python3 utils/generate_ic.py --epoch "${i}"

  SAVE_PATH="data-training/burgers_data_epoch_${i}.npz"

  # Run the monolithic solver and save to save_path
  # The 'None' argument tells the solver to run monolithic (preCICE participant name is none, run without preCICE)
  (
    cd solver-scipy
    ./run.sh --savefile "${SAVE_PATH}"
  )
done

echo "---"
echo "Training data generation complete."
echo "Files are saved in solver-scipy/data-training/"