"""
    AbstractPrescribedStrain

Abstract supertype for macroscopic strain boundary conditions.

See [`MacroscopicStrain`](@ref) for a single load step and
[`MacroscopicStrainSteps`](@ref) for a sequence of load steps.
"""
abstract type AbstractPrescribedStrain end

# Infer spatial dimension from Voigt vector length.
dim_from_length(::Val{3}) = 2
dim_from_length(::Val{6}) = 3
dim_from_length(::Val{N}) where {N} = throw(ArgumentError("Vector must be of length 3 or 6, got $N"))

"""
    MacroscopicStrain{dim, T}(data)
    MacroscopicStrain(data)

Macroscopic strain boundary condition for a single load step.

`data` is a Voigt vector of length 3 (2D) or 6 (3D). When constructed without explicit `dim`,
the spatial dimension is inferred from the vector length.

# Example
```julia
bc = MacroscopicStrain([0.01, 0.0, 0.0, 0.0, 0.0, 0.0])  # uniaxial strain in 3D
bc = MacroscopicStrain{3}([0.01, 0.0, 0.0, 0.0, 0.0, 0.0])  # explicit dim
```
"""
struct MacroscopicStrain{dim, T <: AbstractFloat} <: AbstractPrescribedStrain
    data::Vector{T}

    function MacroscopicStrain{dim, T}(data::Vector{T}) where {dim, T <: AbstractFloat}
        dim in (2, 3) && isa(dim, Integer) || throw(ArgumentError("dim must be 2 or 3"))
        n_voigt = dim == 2 ? 3 : 6
        length(data) == n_voigt || throw(ArgumentError(
            "MacroscopicStrain for dim=$dim requires $n_voigt components, got $(length(data))"
        ))
        new{dim, T}(data)
    end
end

MacroscopicStrain{dim}(data::Vector{T}) where {dim, T <: AbstractFloat} = MacroscopicStrain{dim, T}(data)
MacroscopicStrain(data::Vector{T}) where {T <: AbstractFloat} = MacroscopicStrain{dim_from_length(Val(length(data))), T}(data)

"""
    MacroscopicStrainSteps{dim, T}(steps)

Macroscopic strain boundary condition for a sequence of load steps.

!!! warning
    Not yet supported by any solver. Reserved for future use with nonlinear
    [`HistoryDependent`](@ref) materials.
"""
struct MacroscopicStrainSteps{dim, T <: AbstractFloat} <: AbstractPrescribedStrain
    steps::Vector{Vector{T}}
end

"""
    initialize_strain_field(macro_strain::MacroscopicStrain, microstructure)

Initialize a strain field by filling every voxel with the macroscopic strain `macro_strain`.
"""
initialize_strain_field(
    macro_strain::MacroscopicStrain{dim, T},
    microstructure::InternalMicrostructure{dim, T}
) where {dim, T} = repeat(macro_strain.data, 1, size(microstructure)...)