# 1D Partitioned Burgers' Equation Tutorial

Solving the partitioned 1D Burgers' equation using preCICE with a Neural Network surrogate model.

The monolithic reference solution is computed using a first-order finite volume solver implemented in SciPy.

The partitioned case consists of two participants:
- **Dirichlet**: Solves the left half of the domain.
- **Neumann**: Solves the right half of the domain.

This tutorial includes two versions for the Neumann participant:
1.  A standard finite volume solver (`neumann-scipy`).
2.  A pre-trained neural network surrogate that approximates the solver (`neumann-surrogate`).

## Setup

-  **Create the environment and install dependencies:**

    A script is provided to create and activate the Python virtual environment and install the required packages from `requirements.txt`.

    ```bash
    source activate-env.sh
    ```

## Running the Simulations

- **Run with the SciPy solver for both participants:**

  ```bash
  ./run-partitioned-scipy.sh [seed]
  ```
  Example: `./run-partitioned-scipy.sh 42`

- **Run with the neural network surrogate for the Neumann participant:**

  ```bash
  ./run-partitioned-surrogate.sh [seed]
  ```
  Example: `./run-partitioned-surrogate.sh 42`

Both scripts will execute the partitioned simulation and then compute the monolithic reference solution for comparison.

## Visualizing the Results

The `visualize_partitioned_domain.py` script generates plots comparing the partitioned and monolithic solutions.

- **Generate plots:**

  You can specify which timestep to plot (default is 10).

  ```bash
  python3 visualize_partitioned_domain.py --neumann neumann-surrogate/surrogate.npz [timestep]
  ```
  Example: `python3 visualize_partitioned_domain.py --neumann neumann-surrogate/surrogate.npz 20`

The script will produce the following output files:
- `full_domain_timestep_slice.png`: A plot of the solution u.
- `gradient_timestep_slice.png`: A plot of the gradient du/dx.
- `full_domain_evolution.png`: A plot showing the time evolution over the entire domain.
