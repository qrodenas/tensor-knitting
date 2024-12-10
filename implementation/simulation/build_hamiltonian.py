from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler import CouplingMap

def build_hamiltonian(num_qubits, coupling_map: CouplingMap, decrease = False, initial_qubit_decrease = None, single_gates=['Z'], two_gates=['XX'], anisotropy=1, h=1):
    """
    Builds the Hamiltonian for the quantum system.

    Args:
        num_qubits (int): Total number of qubits.
        coupling_map (CouplingMap): Coupling map defining qubit connectivity.
        decrease (bool): Flag determining if we use coupling decreasing with distance.
        initial_qubit_decrease(int): Initial qubit in which the decreasing is applied. Decreasing will only affect qubits with higher index.
        single_gates (list): List of single-qubit Pauli strings.
        two_gates (list): List of two-qubit Pauli strings.
        anisotropy (float): Anisotropy parameter.
        h (float): Magnetic field strength.

    Returns:
        SparsePauliOp: The constructed Hamiltonian.
    """
    edge_list = coupling_map.get_edges()
    hamlist = []

    for edge in edge_list:
        qubit1, qubit2 = edge
        factor = anisotropy ** (qubit1 - initial_qubit_decrease + 1) if decrease and initial_qubit_decrease is not None and qubit1 >= initial_qubit_decrease else anisotropy
        for gate in two_gates:
            hamlist.append((gate, edge, factor))

    for qubit in coupling_map.physical_qubits:
        for gate in single_gates:
            hamlist.append((gate, [qubit], h))

    hamiltonian = SparsePauliOp.from_sparse_list(hamlist, num_qubits=num_qubits)

    return hamiltonian