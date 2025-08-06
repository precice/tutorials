import os
from bayesvalidrox import PyLinkForwardModel, Input, PCE, ExpDesigns, Engine
import h5py
import joblib
import numpy as np


def create_snapshots() -> None:
    """
    Create snapshots of the DuMuX model in micro-dumux/ using the Micro Manager.
    The snapshots are saved in the specified directory.
    """
    concentration_samples = np.linspace(0.0, 0.5, 50)

    with h5py.File("input_samples.hdf5", "w") as f:
        f.create_dataset("concentration", data=concentration_samples)

    os.subprocess.run(["./create-snapshots.sh", "-s"], check=True)


def read_snapshots(snapshots_dir: str) -> np.array, dict:
    inputs = np.empty((0, 0))
    outputs = {"porosity": [], "k_00": [], "k_01": [], "k_10": [], "k_11": [], "x_values": []}  # Need to be 2D arrays
    return inputs, outputs


def create_surrogate(snapshots_dir: str) -> None:
    # We create a fake model from model.py because we directly provide the
    # input and outputs from the Micro Manager snapshots.
    model = PyLinkForwardModel()
    model.py_file = "model"
    model.name = "micro-dumux-surrogate"

    inputs, outputs = read_snapshots(snapshots_dir)

    inputs = Input()
    inputs.add_marginals("concentration", "unif", )

    exp_design = ExpDesigns()
    exp_design.x = inputs
    exp_design.y = outputs

    # Create the surrogate model
    meta_model = PCE(inputs)
    meta_model.meta_model_type = "aPCE"
    meta_model.pce_reg_method = "FastARD"
    meta_model.pce_degree = 5

    # Train the surrogate model
    engine = Engine(meta_model, model, exp_design)
    engine.train_normal()

    with open(f'{model.name}.pkl', 'wb') as output:
        joblib.dump(engine, output, 2)


def main():
    snapshots_dir = "snapshots"
    if not os.path.exists(snapshots_dir):
        os.makedirs(snapshots_dir)

    os.chdir(snapshots_dir)
    create_snapshots()
    create_surrogate(snapshots_dir)


if __name__ == "__main__":
    main()
