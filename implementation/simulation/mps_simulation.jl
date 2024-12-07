module MPSSimulation

using ITensors

include("gates.jl")
using .Gates
include("../utils/utils.jl")
using .Utils
include("mpo_simulation.jl")
using .MPOSimulation

export mpo_sequence_apply_mps

function mpo_sequence_apply_mps(circ_data, basis_gates, cutoff, maxdim, method)
    all_mpo_sequences = Dict{String, Vector{Any}}()
    all_psis = Dict{String, Vector{Vector{Tuple{MPS, Float64}}}}()
    observables = Dict{String, Vector{MPO}}()
    Zmpo_A = Vector{MPO}()
    Zmpo_B = Vector{MPO}()

    for label in keys(circ_data)
        n_qubits = circ_data[label]["Qubit number"]
        sites_label = create_sites(n_qubits)

        initial_psi = productMPS(sites_label, "0")
        initial_psi_p = Vector{Tuple{MPS, Float64}}()
        push!(initial_psi_p, (initial_psi, 1.0))  # Push the initial state to the MPS list.

        label_mpo_sequences = Vector{Any}()
        all_psis[label] = Vector{Vector{Tuple{MPS, Float64}}}()
        
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

            current_psis = initial_psi_p

            for mpo in mpo_sequence
                if typeof(mpo) != MPO && mpo[1] == "measure"
                    qubit_index = mpo[2]

                    P0_mpo = proj_0_gate(sites_label, qubit_index)
                    P1_mpo = proj_1_gate(sites_label, qubit_index)

                    new_psis = Vector{Tuple{MPS, Float64}}()

                    for (psi, prob) in current_psis
            
                        prob0 = real(inner(psi, apply(P0_mpo, psi)))
                        prob1 = real(inner(psi, apply(P1_mpo, psi)))
                       
                        epsilon = 1e-12  # Define a threshold to approximate zero in measurements
                        if prob0 > epsilon
                            psi0 = apply(P0_mpo, psi)
                            psi0 = normalize(psi0)
                            push!(new_psis, (psi0, prob * prob0))
                        end

                        if prob1 > epsilon
                            psi1 = apply(P1_mpo, psi)
                            psi1 = normalize(psi1)
                            push!(new_psis, (psi1, prob * (-prob1)))
                        end

                    end

                    current_psis = new_psis

                elseif isa(mpo, MPO)
                    new_psis = Vector{Tuple{MPS, Float64}}()

                    for (psi, prob) in current_psis
                        new_psi = apply(mpo, psi; cutoff=cutoff, maxdim=maxdim)
                        push!(new_psis, (new_psi, prob))
                    end

                    current_psis = new_psis
                end
            end

            println("Completed applying MPO sequence for subcircuit $label$index.")

            push!(all_psis[label], current_psis)
        end

        all_mpo_sequences[label] = label_mpo_sequences
    end

    observables["A"] = Zmpo_A
    observables["B"] = Zmpo_B

    return (all_mpo_sequences, all_psis, observables)
end

end 
