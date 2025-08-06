import os
import subprocess
from bayesvalidrox import PyLinkForwardModel, Input, PCE, ExpDesigns, Engine
import h5py
import joblib
import numpy as np


def create_snapshots() -> None:
    """
    Create snapshots of the DuMuX model in micro-dumux/ using the Micro Manager.
    The snapshots are saved in the specified directory.
    """
    # What effect do uniformely distributed concentration samples have on the quality of the surrogate model?
    concentration_samples = np.linspace(0.0, 0.5, 50)

    with h5py.File("input_samples.hdf5", "w") as f:
        f.create_dataset("concentration", data=concentration_samples)

    subprocess.run(["./create-snapshots.sh", "-s"], check=True)


def read_snapshots(snapshots_dir: str) -> tuple:
    """
    Read the snapshots created by the Micro Manager and return the inputs and outputs.
    The inputs are the concentration samples, and the outputs are the porosity and conductivity matrix values.

    Parameters
    ----------
    snapshots_dir : str
        The directory where the snapshots are stored.

    Returns
    -------
    tuple
        A tuple containing the inputs (concentration samples) and outputs (porosity and conductivity values).
    """
    with h5py.File(os.path.join(snapshots_dir, "output.hdf5"), "r") as f:
        concentration_data = f["concentration"]
        porosity_data = f["porosity"]
        k_00_data = f["k_00"]
        k_01_data = f["k_01"]
        k_10_data = f["k_10"]
        k_11_data = f["k_11"]

    concentration = np.swapaxes(np.array([concentration_data]), 0, 1)

    porosity = np.swapaxes(np.array([porosity_data]), 0, 1)
    k_00 = np.swapaxes(np.array([k_00_data]), 0, 1)
    k_01 = np.swapaxes(np.array([k_01_data]), 0, 1)
    k_10 = np.swapaxes(np.array([k_10_data]), 0, 1)
    k_11 = np.swapaxes(np.array([k_11_data]), 0, 1)

    outputs = {"porosity": porosity, "k_00": k_00, "k_01": k_01, "k_10": k_10, "k_11": k_11, "x_values": np.array([0])}
    return concentration, outputs


def split_samples(X, y, n_valid):
    """
    Split the samples and evaluations into training and validation/test data.
    The split is performed randomly.

    Parameters
    ----------
    X : np.ndarray
        Samples, shape (#samples, #parameters)
    y : dict
        Corresponding model evaluations. Expected to match BVR output format.
    n_valid : int
        Number of samples to keep for validation.

    Returns
    -------
    X_train, y_train,
    X_valid, y_valid
    """
    n_samples = X.shape[0]
    if n_valid >= n_samples:
        raise AttributeError('The set number of validation points is invalid.')

    # Random split
    n_train = n_samples - n_valid
    choice = np.random.choice(
        range(n_samples), size=(n_train,), replace=False
    )
    ind = np.zeros(n_samples, dtype=bool)
    ind[choice] = True

    # Split samples
    X_train = X[ind]
    X_valid = X[~ind]

    # Split outputs
    y_train = {}
    y_valid = {}
    for key in y:
        if key != "x_values":
            y_train[key] = y[key][ind]
            y_valid[key] = y[key][~ind]

    return X_train, y_train, X_valid, y_valid


def create_surrogate(snapshots_dir: str) -> None:
    # We create a fake model from model.py because we directly provide the
    # input and outputs from the Micro Manager snapshots.
    model = PyLinkForwardModel()
    model.py_file = "model"
    model.name = "micro-dumux-surrogate"
    model.link_type = "function"
    model.output_names = ["porosity", "k_00", "k_01", "k_10", "k_11"]

    x, y = read_snapshots(snapshots_dir)

    inputs = Input()
    inputs.add_marginals(name="concentration", dist_type="unif", parameters=[0, 0.5])

    exp_design = ExpDesigns(inputs)
    exp_design.x = x
    exp_design.y = y

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
