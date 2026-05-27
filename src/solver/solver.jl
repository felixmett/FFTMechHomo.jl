abstract type AbstractSolver end

abstract type AbstractLinearity end
abstract type Linear <: AbstractLinearity end
abstract type NonLinear <: AbstractLinearity end

struct GridConstants{dim, T <: AbstractFloat}
    n_voigt::Int
    n_cells::T
    zero_idx::CartesianIndex{dim}
end

function GridConstants(microstructure::InternalMicrostructure{dim, T}) where {dim, T}
   n_voigt = dim^2 - dim^(dim-2)
   n_cells = T(prod(size(microstructure.materials)))
   zero_idx = CartesianIndex(ntuple(x -> 1, dim))
   return GridConstants{dim, T}(n_voigt, n_cells, zero_idx)
end


# # Somehow make this aware if linear or not
# struct Solution{S <: AbstractSolver}
#     stress
#     strain
#     stress_avg
#     residuals
# end