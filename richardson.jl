module Richardson
export Mode, modes, αs, rundirname, richardson

struct Mode
    name::AbstractString
    nx::Int
    nd::Int
    nt::Int
end

modes = [Mode("A", 63, 50, 38), Mode("B", 80, 63, 48), Mode("C", 100, 80, 60)]
αs = [0, 4, 8, 10, 12]

function rundirname(run)
    return run[1].name * "alpha" * string(run[2])
end

function richardson_order(values)
    if values[1] .> values[2] .> values[3] .|| values[1] .< values[2] .< values[3]
        return log.((values[1] .- values[2]) ./ (values[2] .- values[3])) ./ log(2)
    else
        return missing
    end
end

# Returns extrapolated value
function richardson(values, n)
    return (2 .^ n .* values[3] .- values[2]) ./ (2 .^ n .- 1)
end

end