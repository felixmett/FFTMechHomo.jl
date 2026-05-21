"""
    AbstractMaterial

Abstract supertype for all material models.

# Voigt Convention
All materials use engineering shear strains (γ = 2ε):

    strain = [ε₁₁, ε₂₂, (ε₃₃,) γ₁₂, (γ₁₃, γ₂₃)]
    stress = [σ₁₁, σ₂₂, (σ₃₃,) σ₁₂, (σ₁₃, σ₂₃)]

# Interface
Every material must implement:
- `compute_stress!(stress, strain, material, i)`
"""
abstract type AbstractMaterial end