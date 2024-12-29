include("richardson.jl")
using .Richardson
using CSV
using DataFrames
using Plots

function getfoam_ts(α, mode)
    basedir = mode * "alpha" * string(α) * "/postProcessing/sampleDict/"
    return readdir(basedir)
end

function foam_cps(α, t, mode)
    file = mode * "alpha" * string(α) * "/postProcessing/sampleDict/" * t * "/CAD.xy"
    foam = DataFrame(CSV.File(file, header=false, skipto=2, delim=" ", ignorerepeated=true))

    xs = foam[:, 1]
    ps = foam[:, 4] / (0.5 * 1.225 * 10 * 10)
    # We compute magnitude
    shear = sum(map(x -> x^2, Matrix(foam[:, 5:10])), dims=2)

    return xs, ps, shear
end

function xfoil_cps(α)
    xfoilfile = "xfoil/alpha" * string(α) * ".csv"
    xfoil = DataFrame(CSV.File(xfoilfile, header=6))

    xs = xfoil[:, "x"]
    ps = xfoil[:, "Cpv"]

    return xs, ps
end

function lastfoam_cps(α, mode)
    return foam_cps(α, getfoam_ts(α, mode)[end], mode)
end

function animatefoam_cps(α, mode)
    ts = getfoam_ts(α, mode)
    xfoilx, xfoil_ps = xfoil_cps(α)
    anim = @animate for t in ts
        (x, p) = foam_cps(α, t, mode)
        pl = scatter(x, p, ylimits=(-5.5, 1.1))
        scatter!(pl, xfoilx, xfoil_ps)
    end
    return anim
end

function plotcylinder(α, mode)
    ts = getfoam_ts(α, mode)
    (xfoilx, xfoil_ps) = xfoil_cps(α)
    (x, p, shear) = foam_cps(α, ts[end], mode)
    pl = scatter(xfoilx, xfoil_ps, label="XFoil", markersize=2, markerstrokewidth=0, xlabel="x(m)", ylabel="cP")
    scatter!(pl, x, p, label="OF", markersize=2, markerstrokewidth=0)
end

function plotcylinder_drag(α, mode)
    ts = getfoam_ts(α, mode)
    (x, p, shear) = foam_cps(α, ts[end], mode)
    scatter(x, shear, label="Shear Magnitude", markersize=2, markerstrokewidth=0, xlabel="x(m)", ylabel="Cτ")
end