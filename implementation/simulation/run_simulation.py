from qiskit_addon_cutting.utils.simulation import ExactSampler
from qiskit_addon_cutting import reconstruct_expectation_values
from qiskit_aer import AerSimulator
import numpy as np

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

def run_mps_simulator(subexperiments):

    mps_simulator = AerSimulator(method='matrix_product_state')
    
    results = {
        label: mps_simulator.run(subexperiment).result().get_counts()
        for label, subexperiment in subexperiments.items()
    }
    return results

def run_mps_simulator2(subexperiments):
    
    mps_simulator = EstimatorV2(backend_options = {'method': 'matrix_product_state'})
    
    results = {
        label: mps_simulator.run(subexperiment).result().get_counts()
        for label, subexperiment in subexperiments.items()
    }
    return results



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