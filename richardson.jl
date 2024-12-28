struct Mode
    name::AbstractString
    nx::Int
    nd::Int
    nt::Int
end

modes = [Mode("A", 60, 40, 40), Mode("B", 85, 56, 56), Mode("C", 120, 80, 80)]
αs = [0, 4, 8, 10, 12]

function rundirname(run)
    return run[1].name * "alpha" * string(run[2])
end