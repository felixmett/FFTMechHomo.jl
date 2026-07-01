using FFTMechHomo
using LinearAlgebra
using Plots

dim = 2
L = 257
r = 64
centered_range(L::Int) = collect(1:L) .- (L÷2 + mod(L,2))

X = centered_range(L)
indicator = [x^2 + y^2 <= r^2 ? 1 : 2 for x in X, y in X]

heatmap(
    indicator;
    color=:jet, colorbar=false,
    xticks=[1, 64, 128, 192, 256],
    yticks=[1, 64, 128, 192, 256],
    xlabel="X", ylabel="Y"
)

mat1 = LinearIsotropicElastic{dim}(72000., 0.22)
mat2 = LinearIsotropicElastic{dim}(48000., 0.3)
microstructure = Microstructure(map(x -> x == 1 ? mat1 : mat2, indicator))

disc = MoulinetSuquetDiscretization(microstructure)
α₀ = 1.5(mat1.μ + mat2.μ)
solver = BasicScheme(α₀, microstructure)

ε0 = 0.01
macro_strains = [
    MacroscopicStrain([ε0, 0, 0]),
    MacroscopicStrain([0, ε0, 0]),
    MacroscopicStrain([0, 0, ε0])
]

ℂ = zeros(3,3)
for (i, macro_strain) in enumerate(macro_strains)
    sol = solve(microstructure, disc, macro_strain, solver)
    @assert sol.converged "Solver has not converged!"

    ℂ[i,:] = sol.stress_avg ./ ε0
end

test_strain = MacroscopicStrain([ε0, ε0, ε0])
test_sol = solve(microstructure, disc, test_strain, solver)
homogenized_stress = ℂ * test_strain.data
isapprox(test_sol.stress_avg, homogenized_stress; rtol=1e-3)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
