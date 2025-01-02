module MPSSimulation

using ITensors, ITensorMPS

include("gates.jl")
using .Gates
include("../utils/utils.jl")
using .Utils
include("mpo_simulation.jl")
using .MPOSimulation

export build_all_mpo_sequences, apply_mpo_sequences, create_mpo, compute_expval, mps_tensor_product, compute_expval_knitting, CoeffExpPair, Results

# Define structures that will be handy for handling the data
struct CoeffExpPair
    coeff::Float64
    exps::Vector{Float64}
end

struct Results
    time::Float16
    expval_true::Float64
    expval_tensors::Float64
    discrepancy::Float64
    maxdim::Int64
end

function build_all_mpo_sequences(circ_data, basis_gates)
    all_mpo_sequences = Dict{String, Vector{Any}}()
    
    for label in keys(circ_data)
        n_qubits = circ_data[label]["Qubit number"]
        sites_label = create_sites(n_qubits)

        label_mpo_sequences = Vector{Any}()

        for subcircuit in circ_data[label]["Subcircuits"]
            operations = subcircuit["Operations"]
            mpo_sequence = build_mpo_sequence(sites_label, operations, basis_gates)
            push!(label_mpo_sequences, mpo_sequence)
        end

        all_mpo_sequences[label] = label_mpo_sequences
    end

    return all_mpo_sequences
end

function apply_mpo_sequences(circ_data, all_mpo_sequences, cutoff, maxdim, method)
    all_psis = Dict{String, Vector{Vector{Tuple{MPS, Float64}}}}()

    for label in keys(circ_data)
        n_qubits = circ_data[label]["Qubit number"]
        sites_label = vcat(siteinds(all_mpo_sequences[label][1][1]; plev = 0)...)

        initial_psi = productMPS(sites_label, "0")
        initial_psi_list = [(initial_psi, 1.0)]
        label_psis = Vector{Vector{Tuple{MPS, Float64}}}()

        # Retrieve all subcircuits and their corresponding MPO sequences
        subcircuits = circ_data[label]["Subcircuits"]
        label_mpo_sequences = all_mpo_sequences[label]


        for (i, subcircuit) in enumerate(subcircuits)
            mpo_sequence = label_mpo_sequences[i]
            current_psis = initial_psi_list

            for mpo in mpo_sequence
                if typeof(mpo) != MPO && mpo[1] == "measure"
                    qubit_index = mpo[2]

                    P0_mpo = proj_0_gate(sites_label, qubit_index)
                    P1_mpo = proj_1_gate(sites_label, qubit_index)

                    new_psis = Vector{Tuple{MPS, Float64}}()

                    for (psi, prob) in current_psis
                        prob0 = real(inner(psi', P0_mpo, psi))
                        prob1 = real(inner(psi', P1_mpo, psi))

                        epsilon = 1e-12
                        if prob0 > epsilon
                            psi0 = apply(P0_mpo, psi; cutoff=cutoff, maxdim=maxdim)
                            psi0 = normalize(psi0)
                            push!(new_psis, (psi0, prob * prob0))
                        end

                        if prob1 > epsilon
                            psi1 = apply(P1_mpo, psi; cutoff=cutoff, maxdim=maxdim)
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

            # println("Completed applying MPO sequence for subcircuit $(label)$(subcircuit["Subcircuit"]).")

            # Store the final states for this subcircuit
            push!(label_psis, current_psis)
        end

        # Store all subcircuits’ final states for this label
        all_psis[label] = label_psis
    end

    return all_psis
end


function create_mpo(op_name::String, site::Int, sites::Vector{Index{Int64}})
    os = OpSum()
    add!(os, op_name, site)
    return MPO(os, sites)
end

function compute_expval(psi::MPS, observable::MPO)
    # Compute ⟨psi| O |psi⟩
    expval = ITensor(1.0)
    for n in 1:length(psi)
        expval *= psi[n]
        expval *= observable[n]
        expval *= prime(dag(psi[n]))
    end
    return real(expval[1])
end

# Faster tensor product
function mps_tensor_product(mps_dict::Dict{String, MPS}, sites::Vector{Index{Int64}})
    total_length = sum(length(mps) for mps in values(mps_dict))
    psi = MPS(sites, "0")
    current_index = 1
    for label in keys(mps_dict)
        mps = mps_dict[label]
        for i in 1:length(mps)
            psi[current_index] = mps[i]
            current_index += 1
        end
    end
    return psi
end

function compute_expval_knitting(QC_data::Vector{CoeffExpPair}, MPS_data::Vector{CoeffExpPair})
    if length(QC_data) != length(MPS_data)
        error("Vectors 'QC_data' and 'MPS_data' must have the same length.")
    end

    total_expval = 0.0

    for i in 1:length(QC_data)
        if QC_data[i].coeff != MPS_data[i].coeff
            error("Mismatch in coefficients at index $i: $(QC_data[i].coeff) vs $(MPS_data[i].coeff)")
        end

        coeff = QC_data[i].coeff
        exps_QC = QC_data[i].exps
        exps_MPS = MPS_data[i].exps

        if length(exps_QC) != length(exps_MPS)
            error("Mismatch in 'exps' length at index $i: $(length(exps_QC)) vs $(length(exps_MPS))")
        end
        exps_product = exps_QC .* exps_MPS
        sum_products = sum(exps_product)
        weighted_sum = coeff * sum_products
        total_expval += weighted_sum
    end

    return total_expval
end


end
