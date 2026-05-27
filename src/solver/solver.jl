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

function residual!(
    stress::AbstractArray{T},
    strain::AbstractArray{T},
    strain_prev::AbstractArray{T},
    stress_avg::AbstractVector{T},
    ref::ReferenceMaterial{dim, T},
    consts::GridConstants{dim, T}
) where {dim, T}
    stress .= strain .- strain_prev
    compute_stress_field!(stress, stress, ref)

    inner_product_stress = T(0)
    for i in CartesianIndices(size(stress)[2:end])
        inner_product_stress += tensor_dot(view(stress, :, i), :stress)
    end
    inner_product_stress /= consts.n_cells
    norm_stress_avg = tensor_dot(stress_avg, :stress)
    
    return sqrt(inner_product_stress / norm_stress_avg)
end

# # Somehow make this aware if linear or not
# struct Solution{S <: AbstractSolver}
#     stress
#     strain
#     stress_avg
#     residuals
# end