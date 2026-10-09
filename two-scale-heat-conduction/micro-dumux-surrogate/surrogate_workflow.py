import os
import subprocess
from bayesvalidrox import PyLinkForwardModel, Input, PCE, ExpDesigns, Engine
import h5py
import joblib
import numpy as np
import matplotlib.pyplot as plt


def create_snapshots() -> None:
    """
    Create snapshots of the DuMuX model in micro-dumux/ using the Micro Manager.
    """
    # What effect do uniformly distributed concentration samples have on the quality of the surrogate model?
    concentration_samples = np.linspace(0.0, 0.9, 500)

    # The name of the HDF5 file should be the same as the parameter_file_name in micro-manager-snapshot-config.json
    with h5py.File("macro-concentration-samples.hdf5", "w") as f:
        f.create_dataset("Concentration", data=concentration_samples)

    # Run the Micro Manager to create snapshot database
    subprocess.run(
        ["mpiexec", "-n", "2", "micro-manager-precice", "--snapshot", "micro-manager-snapshot-config.json"],
        check=True,
    )


def read_snapshots() -> tuple:
    """
    Read the snapshots created by the Micro Manager and return the inputs and outputs.
    The inputs are the concentration samples, and the outputs are the porosity and conductivity matrix values.

    Returns
    -------
    tuple
        A tuple containing the inputs (concentration samples) and outputs (porosity and conductivity values).
    """
    with h5py.File("micro-dumux-snapshots.hdf5", "r") as f:
        concentration_data = f["Concentration"][:]
        porosity_data = f["Porosity"][:]
        k_00_data = f["K00"][:]
        k_11_data = f["K11"][:]

    concentration = np.swapaxes(np.array([concentration_data]), 0, 1)
    porosity = np.swapaxes(np.array([porosity_data]), 0, 1)
    k_00 = np.swapaxes(np.array([k_00_data]), 0, 1)
    k_11 = np.swapaxes(np.array([k_11_data]), 0, 1)

    outputs = {"Porosity": porosity, "K00": k_00, "K11": k_11, "x_values": np.array([0])}
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


def create_surrogate() -> tuple:
    """
    Create a surrogate model from the Micro Manager snapshots.

    Returns
    -------
    tuple
        A tuple containing the validation samples and corresponding model evaluations.
    """
    # We create a fake model from model.py because we directly provide the
    # input and outputs from the Micro Manager snapshots.
    model = PyLinkForwardModel()
    model.py_file = "model"
    model.name = "micro-dumux-surrogate"
    model.link_type = "function"
    model.output_names = ["Porosity", "K00", "K11"]

    x, y = read_snapshots()

    # Split the samples into training and validation sets
    n_valid = 200
    x_train, y_train, x_valid, y_valid = split_samples(x, y, n_valid)

    inputs = Input()
    inputs.add_marginals(name="Concentration", dist_type="unif", parameters=[0, 0.5])

    exp_design = ExpDesigns(inputs)
    exp_design.x = x_train
    exp_design.y = y_train

    # Create the surrogate model
    meta_model = PCE(inputs)
    meta_model.meta_model_type = "aPCE"
    meta_model.pce_reg_method = "FastARD"
    meta_model.pce_deg = 5

    # Train the surrogate model
    engine = Engine(meta_model, model, exp_design)
    engine.train_normal()

    with open(f'{model.name}.pkl', 'wb') as output:
        joblib.dump(engine, output, 2)

    return x_valid, y_valid


def validate_surrogate(x_valid, y_valid, model_name="micro-dumux-surrogate.pkl") -> None:
    """
    Validate the surrogate model using the validation samples and outputs.

    Parameters
    ----------
    x_valid : np.ndarray
        Validation samples.
    y_valid : dict
        Corresponding model evaluations for validation samples.
    model_name : str
        Name of the surrogate model file.
    """
    with open(model_name, 'rb') as input:
        engine = joblib.load(input)

    y_metamod, _ = engine.eval_metamodel(x_valid)

    # engine.plot_adapt(y_valid, y_metamod, y_metamod_std, x_valid)

    # Compare predictions with true values
    plt.figure()
    plt.scatter(y_valid["Porosity"], y_metamod["Porosity"])
    plt.xlabel("True Values")
    plt.ylabel("Predictions")
    plt.title(f"Validation: porosity")
    plt.plot([0.4, 1.1], [0.4, 1.1], "k--")
    plt.xlim(0.4, 1.1)
    plt.ylim(0.4, 1.1)
    plt.savefig("surrogate-validation-porosity.pdf", bbox_inches="tight")

    plt.figure()
    plt.scatter(x_valid, y_valid["Porosity"])
    plt.xlabel("Validation Concentration")
    plt.ylabel("Validation Porosity")
    plt.title("IO: Concentration to Porosity")
    plt.savefig("surrogate-validation-concentration-porosity.pdf", bbox_inches="tight")


def main():
    """
    Create snaphots of the DuMuX model, create a surrogate model from the snapshots, and validate the surrogate model.
    """
    create_snapshots()
    x_valid, y_valid = create_surrogate()
    validate_surrogate(x_valid, y_valid)


if __name__ == "__main__":
    main()
