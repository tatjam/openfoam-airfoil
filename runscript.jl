# run with julia -t 4 runscript.jl (or higher number of threads as desired)
include("richardson.jl")
using .Richardson

function savelines(file, lines)
    open(file, "w") do io
        for line in lines
            write(io, line * "\n")
        end
    end
end

function adjustdir(run)
    dir = rundirname(run)
    # Modify meshing parameters and angle of attack
    lines = readlines(dir * "/meshgen.m")
    lines[7] = "alpha = deg2rad(" * string(run[2]) * ");"
    lines[18] = "Nx = " * string(run[1].nx) * ";"
    lines[19] = "ND = " * string(run[1].nd) * ";"
    lines[20] = "NT = " * string(run[1].nt) * ";"
    savelines(dir * "/meshgen.m", lines)

    # Run octave
    basecmd = Cmd(["octave", "meshgen.m"])
    cmd = Cmd(basecmd, dir=dir)
    Base.run(pipeline(cmd, stdout=devnull, stderr=devnull))
end

function runsim(run)
    dir = rundirname(run)
    basecmd = Cmd(["bash", "./borrar"])
    cmd = Cmd(basecmd, dir=dir)
    Base.run(cmd)
    basecmd = Cmd(["bash", "./ejecutar"])
    cmd = Cmd(basecmd, dir=dir)
    Base.run(cmd)
end

runs = collect(Iterators.product(modes, αs))

Threads.@threads for run in runs
    dir = rundirname(run)
    if isdir(dir)
        #rm(dir, recursive=true)
        continue
    else
        cp("sim", dir)
    end
    adjustdir(run)
    @info "Dir $(dir) is prepared to run, starting..."
    runsim(run)
    @info "Dir $(dir) is done!"
end

# Empezo a las 12:56 (approx) acabó a la 01:21
