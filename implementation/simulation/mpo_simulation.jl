module MPOSimulation

using ITensors

include("gates.jl")
using .Gates
include("../utils/utils.jl")
using .Utils

export build_mpo_sequence, mpo_sequence_apply, apply_measurement

function build_mpo_sequence(sites, operations, basis_gates)
    """
    Builds a list of Matrix Product Operators (MPOs) based on the provided operations.
    """
    mpo_list = []

    for op_dict in operations
        gate_type = op_dict["Name"]
        qubits = op_dict["Qubits"]
        angle = get(op_dict, "Angle", NaN) 

        if gate_type in basis_gates
            if gate_type == "h"
                i = qubits[1]
                push!(mpo_list, h_gate(sites, i))

            elseif gate_type == "rx"
                θ = angle
                i = qubits[1]
                push!(mpo_list, rx_gate(sites, i, θ))

            elseif gate_type == "ry"
                θ = angle
                i = qubits[1]
                push!(mpo_list, ry_gate(sites, i, θ))                

            elseif gate_type == "rz"
                θ = angle
                i = qubits[1]
                push!(mpo_list, rz_gate(sites, i, θ))

            elseif gate_type == "rxx"
                φ = angle
                i, j = qubits[1], qubits[2]
                push!(mpo_list, rxx_gate(sites, i, j, φ))

            elseif gate_type == "rzz"
                φ = angle
                i, j = qubits[1], qubits[2]
                push!(mpo_list, rzz_gate(sites, i, j, φ))

            elseif gate_type == "ryy"
                φ = angle
                i, j = qubits[1], qubits[2]
                push!(mpo_list, ryy_gate(sites, i, j, φ))

            elseif gate_type == "cx"
                control, target = qubits[1], qubits[2]
                push!(mpo_list, cx_gate(sites, control, target))
            end

        elseif gate_type == "measure"
            push!(mpo_list, ("measure", qubits[1])) 

        else
            error("Unsupported gate type: $gate_type. Supported gates are: $(join(basis_gates, ", ")), Measure")
        end
    end

    return mpo_list
end


function mpo_sequence_apply(circ_data, basis_gates, cutoff, maxdim, method)
    all_mpo_sequences = Dict{String, Vector{Any}}()
    all_rhos = Dict{String, Vector{MPO}}()
    observables = Dict{String, Vector{MPO}}()
    Zmpo_A = Vector{MPO}()
    Zmpo_B = Vector{MPO}()

    for label in keys(circ_data)
        n_qubits = circ_data[label]["Qubit number"]
        sites_label = create_sites(n_qubits)

        initial_psi = productMPS(sites_label, "0")
        initial_rho = outer(initial_psi, initial_psi')

        label_mpo_sequences = Vector{Any}()
        label_rhos = Vector{MPO}()

        # Lazy way of defining the observables, change later!
        if label == "A"
            for i in 1:n_qubits
                os = OpSum()
                add!(os, "Z", i) 
                mpo = MPO(os, sites_label)
                push!(Zmpo_A, mpo) 
            end

            for i in 1:n_qubits
                os = OpSum()
                add!(os, "Id", i)  
                mpo = MPO(os, sites_label)
                push!(Zmpo_A, mpo) 
            end
        end

        if label == "B"
            for i in 1:n_qubits
                os = OpSum()
                add!(os, "Id", i) 
                mpo = MPO(os, sites_label)
                push!(Zmpo_B, mpo)  
            end

            for i in 1:n_qubits
                os = OpSum()
                add!(os, "Z", i)  
                mpo = MPO(os, sites_label)
                push!(Zmpo_B, mpo) 
            end
        end

        for subcircuit in circ_data[label]["Subcircuits"]
            operations = subcircuit["Operations"]
            index = subcircuit["Subcircuit"]

            mpo_sequence = build_mpo_sequence(sites_label, operations, basis_gates)
            push!(label_mpo_sequences, mpo_sequence)

            # now we apply the MPO sequence
            rho = initial_rho
            for mpo in mpo_sequence
                if typeof(mpo) != MPO && mpo[1] == "measure"
                    qubit_index = mpo[2]
                    rho = apply_measurement(rho, sites_label, qubit_index, cutoff, maxdim, method)
                elseif isa(mpo, MPO)
                    rho = apply(apply(mpo, rho; cutoff = cutoff, maxdim = maxdim), dag(mpo); cutoff = cutoff, maxdim = maxdim)  # Enhanced dagger by Gian!
                end
            end 
            println("Completed applying MPO sequence for subcircuit $label$index.")
            push!(label_rhos, rho)
        end

        all_mpo_sequences[label] = label_mpo_sequences
        all_rhos[label] = label_rhos

    end
    observables["A"] = Zmpo_A
    observables["B"] = Zmpo_B
    return (all_mpo_sequences, all_rhos, observables)
end

function apply_measurement(rho, sites, index, cutoff, maxdim, method)
    M0 = proj_0_gate(sites, index)
    M1 = proj_1_gate(sites, index)
    return apply(apply(M0, rho; cutoff = cutoff, maxdim = maxdim), dag(M0); cutoff = cutoff, maxdim = maxdim) -
           apply(apply(M1, rho; cutoff = cutoff, maxdim = maxdim), dag(M1); cutoff = cutoff, maxdim = maxdim)
end


end 
