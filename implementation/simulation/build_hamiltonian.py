from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler import CouplingMap

def build_hamiltonian(num_qubits, coupling_map: CouplingMap, decrease = False, single_gates=['Z'], two_gates=['XX'], anisotropy=1, h=1):
    """
    Builds the Hamiltonian for the quantum system.

    Args:
        num_qubits (int): Total number of qubits.
        coupling_map (CouplingMap): Coupling map defining qubit connectivity.
        single_gates (list): List of single-qubit Pauli string names.
        two_gates (list): List of two-qubit Pauli string names.
        anisotropy (float): Anisotropy parameter.
        h (float): Magnetic field strength.

    Returns:
        SparsePauliOp: The constructed Hamiltonian.
    """
    edge_list = coupling_map.get_edges()
    hamlist = []

    for i, edge in enumerate(edge_list):
        factor = anisotropy**(-i) if decrease else anisotropy
        for gate in two_gates:
            hamlist.append((gate, edge, factor))

    for qubit in coupling_map.physical_qubits:
        for gate in single_gates:
            hamlist.append((gate, [qubit], h))


    hamiltonian = SparsePauliOp.from_sparse_list(hamlist, num_qubits=num_qubits)
    
    return hamiltonian