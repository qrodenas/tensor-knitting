### Important Code

The main implementation files and essential functions can be found in the [implementation folder](./implementation).

### Instructions for Running Simulations

1. **Run Qiskit Circuit Simulation**  
   First, open and execute the `qiskit_run.ipynb` notebook. Set your desired parameters within this notebook. This will:
   - Run the Qiskit circuit simulation.
   - Generate a `subcircuits.json` file containing the expectation value, and subcircuits information to recreate the simulation.


2. **Run Tensor Network Simulation and Comparison**  
   After running `qiskit_run.ipynb`, proceed to `julia_part.ipynb`. This notebook will:
   - Use the generated `subcircuits.json` file.
   - Run a tensor network simulation using the Matrix Product Operator (MPO) method.
   - Compare the results with the output from the Qiskit simulation.

