include("richardson.jl")
using .Richardson
using CSV
using DataFrames
using Plots
using ReadableRegex

# Plots the OpenFOAM results alongside XFoil
pl = plot()
xfoil = DataFrame(CSV.File("xfoil/NACA 0012_T1_Re0.200_M0.00_N8.3.csv", header=10))
plot!(pl, xfoil[:, "alpha"], xfoil[:, "CL"], label="CL XFoil")

# Load each value from XFoil

# Returns 12-dimensional force vector (force_P, force_Visc, moment_P, moment_Visc)
function extractforces(forcesstr)
    out = Matrix{Float64}(undef, (1, 12))
    negnumb = maybe(["+", "-"]) * zero_or_more(DIGIT)
    floatreg = capture(negnumb * maybe(".") * negnumb * maybe("e") * negnumb)
    regex = one_or_more("(") * floatreg * " " * floatreg * " " * floatreg * one_or_more(")")
    i = 0
    matches = eachmatch(regex, forcesstr)
    for match in matches
        out[1, i*3+1] = parse(Float64, match.captures[1])
        out[1, i*3+2] = parse(Float64, match.captures[2])
        out[1, i*3+3] = parse(Float64, match.captures[3])
        i += 1
    end
    return out
end

function extractforces_mode(mode, α)
    dir = rundirname((mode, α))
    forces = DataFrame(CSV.File(dir * "/postProcessing/forces/0/forces.dat", header=false, skipto=4))

    return vcat([extractforces(x) for x in forces[:, 2]]...)
end

function extract_t_mode(mode, α)
    dir = rundirname((mode, α))
    forces = DataFrame(CSV.File(dir * "/postProcessing/forces/0/forces.dat", header=false, skipto=4))
    return forces[:, 1]
end

function extractforces_alpha(α)
    return map(x -> extractforces_mode(x, α), modes)
end

function extract_t(α)
    return map(x -> extract_t_mode(x, α), modes)
end

function extractvar(modes, var)
    return map(x -> x[:, var], modes)
end


function plot_convergence(α, var, label)
    ts = extract_t(α)
    a = extractforces_alpha(α)
    v = extractvar(a, var)
    pl = plot()
    scatter!(pl, map(x -> x.nx, modes), map(x -> x[end], v), label=label)
    pl2 = plot()
    for mode in eachindex(modes)
        plot!(pl2, ts[mode], v[mode],
            label="Nx=" * string(modes[mode].nx) * ", ND=" * string(modes[mode].nd) * ", NT=" * string(modes[mode].nx))
    end
    return pl2
end

# The problem is that each mode has used different time-steps 
# What we do is linearly interpolate up to the longest one
function interpolatevar(forces, fromts, tots)
    out = zeros(length(tots), size(forces)[2])
    for i in eachindex(tots)
        t = tots[i]
        left = findlast(x -> x <= t, fromts)
        right = findfirst(x -> x >= t, fromts)
        if isnothing(right)
            @warn "Interpolation failed for t = " * string(t)
            out[i, :] .= forces[end, :]
        elseif left == right
            out[i, :] .= forces[right, :]
        else
            interp = (t - fromts[left]) / (fromts[right] - fromts[left])
            out[i, :] .= forces[left, :] .* (1.0 - interp) .+ forces[right, :] .* interp
        end
    end
    return out
end

function homogeneize(α)
    ts = extract_t(α)
    forces = extractforces_alpha(α)
    maxlen = findmax(map(x -> length(x), ts))[2]
    fixed = similar(forces)
    for mode in eachindex(forces)
        fixed[mode] = interpolatevar(forces[mode], ts[mode], ts[maxlen])
    end
    return fixed
end

function plot_richardsonorder(α, var, label)
    ts = extract_t(α)
    a = extractforces_alpha(α)
    v = extractvar(a, var)
    r = richardson_order.(v)
    pl2 = plot()
    return pl2
end