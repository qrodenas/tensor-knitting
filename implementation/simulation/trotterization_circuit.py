
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.synthesis import LieTrotter
from qiskit import QuantumCircuit

def trotterization_circuit(hamiltonian, trotter_reps, t):

    """
    Generates a quantum circuit for the evolution under a Hamiltonian using LieTrotter approximation for the unitary    

    Args:
        hamiltonian (SparsePauliOp): The previously constructed Hamiltonian.
        trotter_reps (int): number of LieTrotter circuits repeated.
        t (float): length of time of the evolution


    Returns:
        QuantumCircuit: Trotterized evolution quantum circuit
    """
    trotterizator = LieTrotter(reps=trotter_reps, insert_barriers=False)
    U = PauliEvolutionGate(operator=hamiltonian, time=t)
    QC = trotterizator.synthesize(U)

    return QC
