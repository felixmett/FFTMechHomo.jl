abstract type AbstractElastic{dim, T <: AbstractFloat} <: HistoryIndependent{dim, T} end

"""
    LinearIsotropicElastic{dim, T}(E, nu)
    LinearIsotropicElastic{dim}(E, nu)

Linear isotropic elastic material (Hooke's law).

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
struct LinearIsotropicElastic{dim, T} <: AbstractElastic{dim, T}
    E::T
    nu::T
    μ::T
    λ::T

    function LinearIsotropicElastic{dim, T}(E::T, nu::T) where {dim, T <: AbstractFloat}
        dim in (2, 3) && isa(dim, Integer) || throw(ArgumentError("dim must be 2 or 3"))
        E > 0 || throw(ArgumentError("E must be positive"))
        -1 < nu < 0.5 || throw(ArgumentError("nu must be in (-1, 0.5)"))
        μ, λ = lame_constants(E, nu)
        new{dim, T}(E, nu, μ, λ)
    end
end

LinearIsotropicElastic{dim}(E::T, nu::T) where {dim, T <: AbstractFloat} = LinearIsotropicElastic{dim, T}(E, nu)
LinearIsotropicElastic{dim}(E::Real, nu::Real) where dim = LinearIsotropicElastic{dim}(promote(float(E), float(nu))...)

"""
    lame_constants(E::AbstractFloat, nu::AbstractFloat)

Return the Lamé constants `μ` and `λ` derived from the Young's modulus `E` and Poisson ratio `nu`.

# Returns
- `μ`: shear modulus
- `λ`: first Lamé constant

# Example
```julia
μ, λ = lame_constants(210e3, 0.3)
```
"""
function lame_constants(E::AbstractFloat, nu::AbstractFloat)
    μ = E / (2 * (1 + nu))
    λ = E * nu / ((1 + nu) * (1 - 2nu))
    return (μ, λ)
end

"""
    compute_stress!(stress::AbstractArray, strain::AbstractArray, mat::LinearIsotropicElastic, i::CartesianIndex)

Compute the stress response of a [`LinearIsotropicElastic`](@ref) material at index `i`
using Hooke's law: σ = 2με + λ tr(ε) 𝕀.
"""
function compute_stress!(
    stress::AbstractArray{T},
    strain::AbstractArray{T},
    mat::LinearIsotropicElastic{dim, T},
    i::CartesianIndex
) where {dim, T <: AbstractFloat}
    tr_strain = sum(strain[1:dim, i])
    stress[1:dim, i] .= 2mat.μ .* strain[1:dim, i] .+ mat.λ * tr_strain
    # Engineering shear strains: gamma = 2epsilon, so shear stress = μ*gamma (no factor 2 needed)
    stress[dim+1:end, i] .= mat.μ .* strain[dim+1:end, i]
    return
end
