module Utils

using JSON
using ITensors

export h_gate, rz_gate, rx_gate, ry_gate, rxx_gate, ryy_gate, rzz_gate, cx_gate,
       proj_0_gate, proj_1_gate, create_sites, inspect_mpo, parse_subcircuits, mpo_sequence,
       mpo_sequence_apply, apply_measurement


function h_gate(sites, i::Int)
    os = OpSum()
    os += "H", i
    return MPO(os, sites)
end

function rz_gate(sites, i::Int, θ::Float64)
    os = OpSum()
    os += cos(θ / 2.0), "I", i
    os += (-im * sin(θ / 2.0)), "Z", i
    return MPO(os, sites)
end

function rx_gate(sites, i::Int, θ::Float64)
    os = OpSum()
    os += cos(θ / 2.0), "I", i
    os += (-im * sin(θ / 2.0)), "X", i
    return MPO(os, sites)
end

function ry_gate(sites, i::Int, θ::Float64)
    os = OpSum()
    os += cos(θ / 2.0), "I", i
    os += (-im * sin(θ / 2.0)), "Y", i
    return MPO(os, sites)
end

function rxx_gate(sites, i::Int, j::Int, ϕ::Float64)
    os = OpSum()
    os += cos(ϕ / 2.0), "I", i, "I", j
    os += -im * sin(ϕ / 2.0), "X", i, "X", j
    return MPO(os, sites)
end

function ryy_gate(sites, i::Int, j::Int, ϕ::Float64)
    os = OpSum()
    os += cos(ϕ / 2.0), "I", i, "I", j
    os += -im * sin(ϕ / 2.0), "Y", i, "Y", j
    return MPO(os, sites)
end

function rzz_gate(sites, i::Int, j::Int, ϕ::Float64)
    os = OpSum()
    os += cos(ϕ / 2.0), "I", i, "I", j
    os += -im * sin(ϕ / 2.0), "Z", i, "Z", j
    return MPO(os, sites)
end

function cx_gate(sites, control::Int, target::Int)
    os = OpSum()
    os += 0.5, "I", control, "I", target
    os += 0.5, "Z", control, "I", target
    os += 0.5, "I", control, "X", target
    os += -0.5, "Z", control, "X", target
    return MPO(os, sites)
end

function proj_0_gate(sites, i::Int)
    os = OpSum()
    os += 0.5, "I", i
    os += 0.5, "Z", i
    return MPO(os, sites)
end

function proj_1_gate(sites, i::Int)
    os = OpSum()
    os += 0.5, "I", i
    os += -0.5, "Z", i
    return MPO(os, sites)
end


function create_sites(n::Int)
    return [Index(2, "Qubit,i$i") for i in 1:n]
end

function inspect_mpo(mpo::MPO, sites::Vector)
    N = length(sites) 

    println("\nManual Inspection of MPO:")
    for site_idx in 1:N
        println("Site $site_idx Operator:")
        println(mpo[site_idx])
    end
end


function parse_subcircuits(json_filename::String)
    expected_value = nothing
    coefs_list = []
    organized_data = Dict{String, Dict{String, Any}}()

    open(json_filename, "r") do file
        json_data = JSON.parse(file)

        expected_value = get(json_data, "Expected value", nothing)
        coefs_list = get(json_data, "Coefficients", [])

        # Extract "Subcircuits"
        subcircuits = get(json_data, "Subcircuits", [])
        if isempty(subcircuits)
            error("No 'Subcircuits' found in JSON data.")
        end

        # Initialize dictionaries to store subexperiments and qubit numbers
        subexp_dict = Dict{String, Vector{Dict{String, Any}}}()
        qubit_num_dict = Dict{String, Int}()

        for (i, entry) in enumerate(subcircuits)
            # Safely extract keys with default values or handle missing keys
            category = get(entry, "Subexperiment", "Unknown")
            qubit_num = get(entry, "Qubit number", nothing)
            qubit_range = get(entry, "Qubit range", [1, 1])
            index = get(entry, "Subcircuit", i - 1)  # Fallback to enumeration if missing
            operations_raw = get(entry, "Operations", [])

            # Validate essential keys
            if qubit_num === nothing
                error("Entry at index $i is missing the 'Qubit number' key.")
            end

            if length(qubit_range) != 2
                error("Entry at index $i has an invalid 'Qubit range'. It should be a two-element array.")
            end

            offset = qubit_range[1] - 1

            if haskey(qubit_num_dict, category)
                if qubit_num_dict[category] != qubit_num
                    error("Inconsistent qubit numbers in category '$category'. Expected $(qubit_num_dict[category]), found $qubit_num.")
                end
            else
                qubit_num_dict[category] = qubit_num
            end

            processed_operations = Vector{Dict{String, Any}}()
            for (j, op) in enumerate(operations_raw)
                name = get(op, "Name", "unknown_operation")
                qubits = get(op, "Qubits", [])
                angle_array = get(op, "Angle", [])

                if isempty(qubits)
                    error("Operation at subcircuit index $index, operation $j is missing the 'Qubits' key.")
                end

                adjusted_qubits = [q - offset for q in qubits]

                angle = (length(angle_array) == 1) ? angle_array[1] : NaN

                processed_op = Dict{String, Any}(
                    "Name" => name,
                    "Angle" => angle,
                    "Qubits" => adjusted_qubits
                )

                push!(processed_operations, processed_op)
            end

            subcircuit = Dict{String, Any}(
                "Subcircuit" => index + 1,  # Adjusting index to start at 1
                "Operations" => processed_operations
            )

            if haskey(subexp_dict, category)
                push!(subexp_dict[category], subcircuit)
            else
                subexp_dict[category] = [subcircuit]
            end
        end

        sorted_categories = sort(collect(keys(subexp_dict)))

        for category in sorted_categories
            sorted_subcircuits = sort(subexp_dict[category], by = sc -> sc["Subcircuit"])
            organized_data[category] = Dict{String, Any}(
                "Qubit number" => qubit_num_dict[category],
                "Subcircuits" => sorted_subcircuits
            )
        end
    end  
    return organized_data, expected_value, coefs_list
end

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


function apply_measurement(rho, sites, index, cutoff, maxdim, method)
    M0 = proj_0_gate(sites, index)
    M1 = proj_1_gate(sites, index)
    return apply(apply(M0, rho; cutoff = cutoff, maxdim = maxdim), dag(M0); cutoff = cutoff, maxdim = maxdim) - apply(apply(M1, rho; cutoff = cutoff, maxdim = maxdim), dag(M1); cutoff = cutoff, maxdim = maxdim)
end


# Combined function to build MPO sequences and apply them to rho
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

        # os = OpSum()
        # for i in 1:n_qubits
        #     os += "Z", i
        # end
        # Zmpo = MPO(os, sites_label)


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
                    rho = apply(apply(mpo, rho; cutoff = cutoff, maxdim = maxdim), swapprime(dag(mpo), 0 => 1); cutoff = cutoff, maxdim = maxdim)  # Enhanced dagger by Gian!
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

end
