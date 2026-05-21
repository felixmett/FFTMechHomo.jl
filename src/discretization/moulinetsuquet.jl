struct MoulinetSuquetDiscretization{dim, T <: AbstractFloat} <: AbstractDiscreteGreenOperator
    grid_size::NTuple{dim, Int}
    ξ::NTuple{dim, Vector{T}}
end

function MoulinetSuquetDiscretization(microstructure::MicroStructure{dim, T}) where {dim, T <: AbstractFloat}
    dim in (2, 3) && isa(dim, Integer) || throw(ArgumentError("dim must be 2 or 3"))
    grid_size = size(microstructure)
    ξ = ntuple(i -> i < dim ? fftfreq(grid_size[i], T(1)) : rfftfreq(grid_size[i], T(1)), dim)
    MoulinetSuquetDiscretization{dim, T}(grid_size, ξ)
end

Base.ndims(::MoulinetSuquetDiscretization{dim}) where {dim} = dim
Base.eltype(::MoulinetSuquetDiscretization{dim, T}) where {dim, T <: AbstractFloat} = T

"""
    Γ⁰!(field, ref, disc::MoulinetSuquetDiscretization)

Moulinec-Suquet discretization of the Green operator.
"""
function Γ⁰!(
    field::AbstractArray{Complex{T}}, 
    ref::ReferenceMaterial{dim,T}, 
    disc::MoulinetSuquetDiscretization{dim,T}
) where {dim, T <: AbstractFloat}
    field_tensor = MMatrix{dim, dim, Complex{T}}(undef)
    voigt_idx = voigt_index(Val(dim))

    fft_grid_size = ntuple(i -> length(disc.ξ[i]), dim)
    for idx in CartesianIndices(fft_grid_size)
        if is_zero_or_nyquist(disc, idx)
            field[:, idx] .= zero(Complex{T})
            continue
        end

        k = SVector(ntuple(i -> disc.ξ[i][idx[i]], dim))
        eta = normalize(k)
        
        for i in 1:dim, j in 1:dim
            field_tensor[i,j] = field[voigt_idx[i,j], idx]
        end
        
        f = field_tensor * eta
        s = dot(f, eta)
        u = (1 / 0.5ref.α₀) * (-f + T(0.5) * s * eta)

        field[:, idx] .= zero(Complex{T})
        for i in 1:dim, j in 1:dim
            field[voigt_idx[i,j], idx] += u[i] * eta[j]
        end
    end
end