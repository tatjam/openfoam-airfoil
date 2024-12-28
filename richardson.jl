module Richardson
export Mode, modes, αs, rundirname, richardson

struct Mode
    name::AbstractString
    nx::Int
    nd::Int
    nt::Int
end

#modes = [Mode("A", 50, 40, 30), Mode("B", 71, 57, 42), Mode("C", 100, 80, 60)]
#αs = [0, 4, 8, 10, 12]

modes = [Mode("A", 50, 40, 30), Mode("B", 71, 57, 42), Mode("C", 100, 80, 60)]
αs = [8]

function rundirname(run)
    return run[1].name * "alpha" * string(run[2])
end

# Returns order and extrapolated value
function richardson(values)
    n = log((values[1] - values[2]) / (values[2] - values[3])) / log(2)
    return n, (2^n * values[3] - values[2]) / (2^n - 1)
end

end