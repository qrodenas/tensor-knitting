
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.synthesis import LieTrotter

def trotterization_circuit(hamiltonian, trotter_reps, dt):

    """
    Generates a quantum circuit for the evolution under a Hamiltonian using LieTrotter approximation for the unitary    

    Args:
        hamiltonian (SparsePauliOp): The previously constructed Hamiltonian.
        trotter_reps (int): number of LieTrotter circuits repeated.
        dt (float): length of time of the evolution


    Returns:
        QuantumCircuit: Trotterized evolution quantum circuit
    """
    trotterizator = LieTrotter(reps=trotter_reps, insert_barriers=False)
    U = PauliEvolutionGate(operator=hamiltonian, time=dt)
    QC = trotterizator.synthesize(U)

    return QC