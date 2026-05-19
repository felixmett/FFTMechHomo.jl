abstract type AbstractDiscreteGreenOperator end

function is_first_or_nyquist(idx::CartesianIndex, disc::AbstractDiscreteGreenOperator)
    all(Tuple(idx) .== 1) && return true
    grid_size = disc.grid_size
    for (dim_idx, i) in enumerate(Tuple(idx))
        n = grid_size[dim_idx]
        n % 2 == 0 && i == n ÷ 2 && return true
    end
    return false
end