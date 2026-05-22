"""
    ReferenceMaterial{dim, T}(α₀)
    ReferenceMaterial{dim}(α₀)

Linear isotropic elastic reference material (Hooke's law).
The reference material tensor is assumed to have the form ℂ=α₀𝕀,
i.e. Lamé parameters μ=0.5α₀ and λ=0.

See [`AbstractMaterial`](@ref) for Voigt convention.

# Type Parameters
- `dim`: spatial dimension (2 or 3)
- `T`: numeric type, must be `<: AbstractFloat`

# Arguments
- `α₀`: material stiffness parameter, must be positive. Related to the shear modulus by α₀ = 2μ.

# Example
```julia
ref = ReferenceMaterial{3}(1.0)
```
"""
struct ReferenceMaterial{dim, T <: AbstractFloat} <: AbstractElastic
    α₀::T

    function ReferenceMaterial{dim, T}(α₀::T) where {dim, T <: AbstractFloat}
        dim in (2, 3) && isa(dim, Integer) || throw(ArgumentError("dim must be 2 or 3"))
        α₀ > 0 || throw(ArgumentError("α₀ must be positive"))
        new{dim, T}(α₀)
    end
end

ReferenceMaterial{dim}(α₀::T) where {dim, T <: AbstractFloat} = ReferenceMaterial{dim, T}(α₀)
ReferenceMaterial{dim}(α₀::Real) where dim = ReferenceMaterial{dim}(float(α₀))

"""
    compute_stress!(stress::AbstractArray, strain::AbstractArray, mat::ReferenceMaterial, i::CartesianIndex)

Compute the stress response of a [`ReferenceMaterial`](@ref), which corresponds
to a [`LinearIsotropicElastic`](@ref) material with μ=0.5α₀ and λ=0.

See [`AbstractMaterial`](@ref) for Voigt convention.
"""
function compute_stress!(
    stress::AbstractArray{T},
    strain::AbstractArray{T},
    mat::ReferenceMaterial{dim, T},
    i::CartesianIndex
) where {dim, T <: AbstractFloat}
    stress[1:dim, i] .= mat.α₀ .* strain[1:dim, i]
    # Engineering shear strains from Voigt notation: γ = 2ε, so τ = 0.5α₀γ (μ instead of 2μ)
    stress[dim+1:end, i] .= 0.5mat.α₀ .* strain[dim+1:end, i]
    return
end

"""
    subtract_stress!(stress::AbstractArray, strain::AbstractArray, mat::ReferenceMaterial)

Subtract the stress response of a [`ReferenceMaterial`](@ref) from `stress` in-place at every field index.
Unlike [`compute_stress!`](@ref), this function operates over the full field rather than at a single index.

See [`AbstractMaterial`](@ref) for Voigt convention.
"""
function subtract_stress!(
    stress::AbstractArray{T},
    strain::AbstractArray{T},
    mat::ReferenceMaterial{dim, T}
) where {dim, T <: AbstractFloat}
    for i in CartesianIndices(size(stress)[2:end])
        stress[1:dim, i] .-= mat.α₀ .* strain[1:dim, i]
        stress[dim+1:end, i] .-= 0.5mat.α₀ .* strain[dim+1:end, i]
    end
    return
end