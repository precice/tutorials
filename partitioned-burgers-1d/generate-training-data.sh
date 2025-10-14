#!/usr/bin/env bash
set -e -u

mkdir -p solver-scipy-fvolumes/data-training

# Number of training runs to generate
NUM_RUNS=200

echo "Generating ${NUM_RUNS} training data samples..."

for i in $(seq 0 $((${NUM_RUNS}-1)))
do
  echo "--- Generating epoch ${i} ---"
  
  # Generate IC
  python3 generate_ic.py --epoch ${i}

  SAVE_PATH="data-training/burgers_data_epoch_${i}.npz"

  # Run the monolithic solver and save to save_path
  # The 'None' argument tells the solver to run monolithic (preCICE participant name is none, run without preCICE)
  cd solver-scipy-fvolumes
  python3 solver.py None --savefile "${SAVE_PATH}"
  cd ..
done

echo "---"
echo "Training data generation complete."
echo "Files are saved in solver-scipy-fvolumes/data-training/"
