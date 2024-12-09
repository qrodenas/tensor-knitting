from qiskit_addon_cutting.utils.simulation import ExactSampler
from qiskit_addon_cutting import reconstruct_expectation_values
from qiskit_aer import AerSimulator
import numpy as np
from qiskit_aer.primitives import EstimatorV2
from qiskit.primitives import StatevectorEstimator

def run_exact_sampler(subexperiments):
    """
    Executes ExactSampler (from qiskit knitting addon) on the sub-experiments.

    Args:
        subexperiments (dict): Dictionary of sub-experiments to run.

    Returns:
        dict: Results from each sub-experiment.
    """
    exact_sampler = ExactSampler()
    results = {
        label: exact_sampler.run(subexperiment).result()
        for label, subexperiment in subexperiments.items()
    }
    return results


def run_mps_simulator_knitted(subexperiments, observables, coeffs, shots):
    mps_simulator = EstimatorV2(options={'backend_options': {'method': 'matrix_product_state'}, 'run_options': {'shots': shots}})
    results = {label: [] for label in subexperiments.keys()}
    for label, circuits in subexperiments.items():
        for circuit in circuits:
            job = mps_simulator.run([(circuit, observables[label])])
            exp = job.result()[0].data.evs
            results[label].append(exp)

    expval = 0.0
    for i in range(len(coeffs)):
        product = coeffs[i][0] 
        for label in results:
            product *= results[label][i] 
        expval += np.sum(product)
    return expval


def run_mps_simulator_full(circuit, observables, shots):
    mps_simulator = EstimatorV2(options={'backend_options': {'method': 'matrix_product_state'}, 'run_options': {'shots': shots}})
    
    job = mps_simulator.run([(circuit, observables)])
    expval = job.result()[0].data.evs
    return expval

def run_statevector_simulator_full(circuit, observables):
    estimator = StatevectorEstimator()

    job = estimator.run([(circuit, observables)])
    expval = job.result()[0].data.evs
    return expval


def reconstruct_expectation(results, coefficients, subobservables, z_observables):
    """
    Reconstructs the expectation value from the results

    Args:
        results (dict): Results from each sub-experiment.
        coefficients (list): Coefficients associated with the sub-experiments (and their sub-circuits).
        subobsevables (SparsePauliOp): observables generated that we want to measure.
        z_observables (SparsePauliOp): overall observables we want to measure in the circuit, we want to extract their coefficients.

    Returns:
        float: The reconstructed expectation value.
    """
    
    reconstructed_expval_terms = reconstruct_expectation_values(results, coefficients, subobservables)
    reconstructed_expval = np.dot(reconstructed_expval_terms, z_observables.coeffs)
    return reconstructed_expval