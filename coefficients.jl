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
    out = Array{Float64}(undef, 12)
    negnumb = maybe(["+", "-"]) * zero_or_more(DIGIT)
    floatreg = capture(negnumb * maybe(".") * negnumb * maybe("e") * negnumb)
    regex = one_or_more("(") * floatreg * " " * floatreg * " " * floatreg * one_or_more(")")
    i = 0
    for match in eachmatch(regex, forcesstr)
        out[i*3+1] = parse(Float64, match.captures[1])
        out[i*3+2] = parse(Float64, match.captures[2])
        out[i*3+3] = parse(Float64, match.captures[3])
        i += 1
    end
    return out
end


function richardsonCLα(α)
    # Load each mode 
    for mode in modes
        dir = rundirname((mode, α))
        forces = DataFrame(CSV.File(dir * "/postProcessing/forces/0/forces.dat", header=false, skipto=4))

        display(forces[:, 2])
    end
end

CLs = map(x -> richardsonCLα(x), αs)
