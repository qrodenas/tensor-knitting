module Utils

using JSON, ITensors

export parse_subcircuits, create_sites, inspect_mpo


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



function create_sites(n::Int)
    return [Index(2, "Qubit,i$i") for i in 1:n]
end

function inspect_mpo(mpo::MPO, sites::Vector{Index})
    N = length(sites)

    println("\nManual Inspection of MPO:")
    for site_idx in 1:N
        println("Site $site_idx Operator:")
        println(mpo[site_idx])
    end
end


end 
