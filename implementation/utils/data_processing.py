from qiskit import transpile
import numpy as np

def circuit_data_dict(subexperiments, coefs, expval, basis_gates):
    """
    Processes sub-experiments to extract data and places it in an organized dictionary.

    Args:
        subexperiments (dict): Dictionary of sub-experiments.
        coefs (list): List of coefficients for reconstruction of the expectation value  
        expval (complex): Expectation value of the full circuit  
        basis_gates (list): List of basis gates for transpilation.


    Returns:
        dict: Dictionary containing the expected value of the experiment, the coefficients for reconstruction of the expval and the sub-experiments data.
    """
    readable_data = []
    for category in sorted(subexperiments.keys()):
        for index, circuit in enumerate(subexperiments[category]):
            transpiled = transpile(circuit, basis_gates=basis_gates) # TODO add optimization when transpiling?
            operations = []
            qubit_range = [transpiled.qubits[0]._index, transpiled.qubits[-1]._index]
            
            for instr in transpiled.data:
                name = instr.operation.name
                angle = instr.operation.params
                qubits = [q._index for q in instr.qubits]
                
                if name == 'measure' and (
                    not instr.clbits or
                    instr.clbits[0]._register.name != 'qpd_measurements'
                ):
                    continue
                    
                operations.append({
                    "Name": name,
                    "Angle": angle,
                    "Qubits": qubits
                })
            
            readable_data.append({
                "Subexperiment": category,
                "Subcircuit": index,
                "Qubit number": transpiled.num_qubits,
                "Qubit range": qubit_range,
                "Operations": operations
            })

    coefficients_list = [coefs[i][0] for i in range(len(coefs))]
    output_dict = {
        "Expected value": np.real(expval),
        "Coefficients": coefficients_list,
        "Subcircuits": readable_data
    }
    return output_dict

def save_to_json(output_dict, filename="subcircuits.json"):
    """
    Saves the output dictionary to a JSON file.

    Args:
        output_dict (dict): The data to be saved.
        filename (str): The name of the JSON file.
    """
    import json
    with open(filename, "w") as json_file:
        json.dump(output_dict, json_file, indent=4)
    print(f"Subcircuits data has been successfully saved to '{filename}'.")
