"""
    AbstractMaterial

Abstract supertype for all material models.

# Voigt Convention
All materials use engineering shear strains (γ = 2ε):

    strain = [ε₁₁, ε₂₂, (ε₃₃, γ₂₃, γ₁₃,) γ₁₂]
    stress = [σ₁₁, σ₂₂, (σ₃₃, σ₂₃, σ₁₃,) σ₁₂]

Subtypes are categorized by memory behaviour:
- [`HistoryIndependent`](@ref): stress depends only on current strain
- [`HistoryDependent`](@ref): stress depends on strain history and internal variables
"""
abstract type AbstractMaterial end

"""
    HistoryIndependent <: AbstractMaterial

Abstract supertype for materials whose stress response depends only on the
current strain, with no internal variables. Example: linear elasticity
"""
abstract type HistoryIndependent <: AbstractMaterial end

"""
    HistoryDependent <: AbstractMaterial

Abstract supertype for materials with internal variables that evolve over the
deformation history. Examples: plasticity, damage, viscoelasticity.
"""
abstract type HistoryDependent <: AbstractMaterial end

"""
    compute_stress!(stress, strain, material::AbstractMaterial, i::CartesianIndex)

Compute the stress response of `material` at grid index `i` in-place.

# Interface
Subtypes of [`AbstractMaterial`](@ref) must implement this method. The stress
and strain arrays follow the Voigt convention defined in [`AbstractMaterial`](@ref).
"""
function compute_stress! end