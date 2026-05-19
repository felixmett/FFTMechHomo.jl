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

"""
    voigt_index(::Val{2})
    voigt_index(::Val{3})

Return the Voigt index mapping for a symmetric matrix of dimension `d`.

The returned nested tuple `vidx` satisfies `vidx[i][j] == v` where `v` is the
index into the Voigt vector corresponding to the `(i, j)` entry of a symmetric
matrix. See [`AbstractMaterial`](@ref) for the Voigt convention used throughout.

# Returns
- `NTuple{d, NTuple{d, Int}}`: Nested tuple representing the `d×d` index matrix

# Example
```julia
vidx = voigt_index(Val(2))
vidx[1][2]  # → 3, i.e. the (1,2) matrix entry is stored at Voigt index 3
```
"""
voigt_index(::Val{2}) = ((1, 3), (3, 2))
voigt_index(::Val{3}) = ((1, 6, 5), (6, 2, 4), (5, 4, 3))