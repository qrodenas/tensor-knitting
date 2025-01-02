module Gates

using ITensors, ITensorMPS

export h_gate, rz_gate, rx_gate, ry_gate, rxx_gate, ryy_gate, rzz_gate, cx_gate,
       proj_0_gate, proj_1_gate

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

end 
