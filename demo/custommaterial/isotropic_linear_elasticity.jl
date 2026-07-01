using FFTMechHomo

struct CustomElasticity{dim, T <: AbstractFloat} <: HistoryIndependent{dim, T}
    E::T
    nu::T
    μ::T
    λ::T
    function CustomElasticity{dim, T}(E::T, nu::T) where {dim, T <: AbstractFloat}
        μ = E / (2 * (1 + nu))
        λ = E * nu / ((1 + nu) * (1 - 2nu))
        new{dim, T}(E, nu, μ, λ)
    end
end

function FFTMechHomo.compute_stress!(
    stress::AbstractArray{T},
    strain::AbstractArray{T},
    mat::CustomElasticity{dim, T},
    i::CartesianIndex
) where {dim, T <: AbstractFloat}
    tr_strain = sum(strain[1:dim, i])
    stress[1:dim, i] .= 2mat.μ .* strain[1:dim, i] .+ mat.λ * tr_strain
    stress[dim+1:end, i] .= mat.μ .* strain[dim+1:end, i]
    return
end

simple_microstructure(mat1, mat2) = Microstructure([i==1 ? mat1 : mat2 for i in 1:3, j in 1:3, k in 1:3])

dim = 3
mat_soft_custom = CustomElasticity{dim, Float64}(2100., 0.3)
mat_hard_custom = CustomElasticity{dim, Float64}(72000., 0.22)
microstructure = simple_microstructure(mat_soft_custom, mat_hard_custom)

disc = MoulinetSuquetDiscretization(microstructure)
macro_strain = MacroscopicStrain([0.01, 0, 0, 0, 0, 0])
solver = BasicScheme(mat_soft_custom.μ + mat_hard_custom.μ, microstructure)
sol = solve(microstructure, disc, macro_strain, solver)

mat_soft_builtin = LinearIsotropicElastic{dim, Float64}(2100., 0.3)
mat_hard_builtin = LinearIsotropicElastic{dim, Float64}(72000., 0.22)
validation_microstructure = simple_microstructure(mat_soft_builtin, mat_hard_builtin)

disc = MoulinetSuquetDiscretization(validation_microstructure)
solver = BasicScheme(mat_soft_builtin.μ + mat_hard_builtin.μ, validation_microstructure)
validation_sol = solve(validation_microstructure, disc, macro_strain, solver)

@assert sol.stress_avg ≈ validation_sol.stress_avg

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
