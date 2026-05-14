abstract type AbstractElastic <: AbstractMaterial end

"""
    LinearIsotropicElastic{dim, T}(E, nu)
    LinearIsotropicElastic{dim}(E, nu)

Linear isotropic elastic material (Hooke's law).

See [`AbstractMaterial`](@ref) for Voigt convention.

# Type Parameters
- `dim`: spatial dimension (2 or 3)
- `T`: numeric type, must be `<: AbstractFloat`

# Arguments
- `E::T`: Young's modulus, must be positive
- `nu::T`: Poisson's ratio, must be in (-1, 0.5)

# Example
```julia
mat = LinearIsotropicElastic{3}(210e3, 0.3) # Float64 by default
mat = LinearIsotropicElastic{3, Float32}(210e3, 0.3) # explicit Float32)
```
"""
struct LinearIsotropicElastic{dim, T <: AbstractFloat} <: AbstractElastic
    E::T
    nu::T

    function LinearIsotropicElastic{dim, T}(E::T, nu::T) where {dim, T <: AbstractFloat}
        dim in (2, 3) && isa(dim, Integer) || throw(ArgumentError("dim must be 2 or 3"))
        E > 0 || throw(ArgumentError("E must be positive"))
        -1 < nu < 0.5 || throw(ArgumentError("nu must be in (-1, 0.5)"))
        new{dim, T}(E, nu)
    end
end

LinearIsotropicElastic{dim}(E::T, nu::T) where {dim, T <: AbstractFloat} = LinearIsotropicElastic{dim, T}(E, nu)
LinearIsotropicElastic{dim}(E::Real, nu::Real) where dim = LinearIsotropicElastic{dim}(promote(float(E), float(nu))...)

"""
    compute_stress!(stress, strain, material::LinearIsotropicElastic, i)

See [`AbstractMaterial`](@ref) for Voigt convention.
"""
function compute_stress!(
    stress::AbstractArray{T},
    strain::AbstractArray{T},
    material::LinearIsotropicElastic{dim, T},
    i::CartesianIndex
) where {dim, T <: AbstractFloat}
    tr_strain = sum(strain[1:dim, i])
    mu = material.E / (2 * (1 + material.nu))
    lambda = material.E * material.nu / ((1 + material.nu) * (1 - 2material.nu))

    stress[1:dim, i] .= 2mu .* strain[1:dim, i] .+ lambda * tr_strain
    # Engineering shear strains: gamma = 2epsilon, so shear stress = mu*gamma (no factor 2 needed)
    stress[dim+1:end, i] .= mu .* strain[dim+1:end, i]
    return
end
