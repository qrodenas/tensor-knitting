from qiskit_addon_cutting import partition_problem, generate_cutting_experiments

def cut_circuit(circuit, partition_labels, observables, num_samples):
    """
    Partitions the given quantum circuit based on the partition labels and generates sub-experiments.

    Args:
        circuit (QuantumCircuit): The quantum circuit to partition.
        partition_labels (str): Labels indicating how to partition the circuit (for 4 qubits, for example AABB).
        observables (SparsePauliOp): Observables for which expectation values are calculated.
        num_samples (int or float): Number of samples for each sub-experiment.

    Returns:
        tuple: Sub-experiments, coefficients, and sub-observables.
    """
    partitioned_problem = partition_problem(
        circuit=circuit,
        partition_labels=partition_labels,
        observables=observables.paulis
    )
    
    subcircuits = partitioned_problem.subcircuits
    subobservables = partitioned_problem.subobservables
    
    subexperiments, coefficients = generate_cutting_experiments(
        circuits=subcircuits,
        observables=subobservables,
        num_samples=num_samples
    )
    return subexperiments, coefficients, subobservables
