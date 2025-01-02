module ReconstructExpectation

using ITensors, ITensorMPS, IterTools

export reconstruct_expectation

function gather_sites(all_psis::Dict{String, Vector{Vector{Tuple{MPS, Float64}}}})
    labels = keys(all_psis)
    sites_per_label = Dict{String, Vector{Index{Int64}}}()
    # Collect site indices for each label
    for label in labels
        sites_per_label[label] = siteinds(all_psis[label][1][1][1])
    end
    # Generate union of sites
    sites_union = vcat(union(values(sites_per_label))...)
    return sites_union
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


function reconstruct_expectation(
    all_psis::Dict{String, Vector{Vector{Tuple{MPS, Float64}}}}, 
    coefficients::Vector{Any}
)

    # We reconstruct the expectation value (with MPS)
    # Tr[Oρ'] = ∑ₓⁿ Tr[Oᴀₓ ⊗ Oᵦₓ ρ'] = ∑ₓⁿ ∑ₖ₌₁ᴷ ∑ⱼₐ₌₁ᵐₐ ∑ⱼᵦ₌₁ᵐᵦ αₖ pᴀⱼₐ pᵦⱼᵦ ⟨Oᴀₓ⟩ₖ,ⱼₐ ⟨Oᵦₓ⟩ₖ,ⱼᵦ

    labels = collect(keys(all_psis))    
    sites_union = gather_sites(all_psis)
    # Create the observable 
    obs = sum([create_mpo("Z", i, sites_union) for i in 1:length(sites_union)]) # ZII + IZI + IIZ
    # obs = [MPO(sites_union, "Z")] # ZZZ

    expval_tensors = 0.0
    for (k, coeff) in enumerate(coefficients)
        temp = 0.0
        combinations = IterTools.product([all_psis[label][k] for label in labels]...)
        for combination in combinations
            mps_dict = Dict{String, MPS}()
            prob = 1.0
            for (idx, label) in enumerate(labels)
                psi_label, p_label = combination[idx]
                mps_dict[label] = psi_label
                prob *= p_label
            end
            psi = mps_tensor_product(mps_dict, sites_union)
            temp += prob * compute_expval(psi, obs)
        end
        expval_tensors += coeff * temp
    end
    return expval_tensors
end

end 
