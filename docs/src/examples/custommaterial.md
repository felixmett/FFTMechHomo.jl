# User-defined History Independent Material

In this example, we define a custom history-independent linear elastic material and integrate it into the `FFTMechHomo.jl` framework. We validate our implementation by comparing it against the built-in isotropic linear elastic material model.
## Isotropic Linear Elasticity

Isotropic linear elasticity is fully characterized by two parameters, which we choose as the Lamé parameters ``\mu`` and ``\lambda``.

The stress–strain relation is given by:

```math
\sigma_{ij} = 2\mu\varepsilon_{ij} + \lambda\varepsilon_{kk}\delta_{ij},
```

where ``\delta_{ij}`` is the Kronecker-delta.

To integrate a custom material into `FFTMechHomo.jl`, we define a new type that implements the required stress update function. Note that we construct the material instance using the Young's modulus ``E`` and the Poisson ratio ``\nu``, that internally get transformed to the Lamé Parameters. We do this to match the built-in material model's constructor for easier comparison.

````@example custommaterial
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
nothing # hide
````

For internal consistency checks, we need to specify a dimension `dim` and a data type `T` of the parameters.

To integrate the new material into `FFTMechHomo.jl`, we implement the required stress update function. Due to the internal use of Voigt-notation, the factor ``2`` for shear stresses is already adsorbed in the shear strains, therefore we only need ``\mu`` for these components.

This function specifies how the material computes stress from a given strain and is called by the solver for each point in the microstructure.

````@example custommaterial
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
nothing # hide
````

With that our material can be used inside of the solver of `FFTMechHomo.jl`. We now test our implementation against the preexisting isotropic linear elasticity model. For that, we first define a simple microstructure.

````@example custommaterial
simple_microstructure(mat1, mat2) = Microstructure([i==1 ? mat1 : mat2 for i in 1:3, j in 1:3, k in 1:3])
nothing # hide
````

Then we create a microstructure with instances of our new material.

````@example custommaterial
dim = 3
mat_soft_custom = CustomElasticity{dim, Float64}(2100., 0.3)
mat_hard_custom = CustomElasticity{dim, Float64}(72000., 0.22)
microstructure = simple_microstructure(mat_soft_custom, mat_hard_custom)
nothing # hide
````

To compute the strain response of the heterogeneous material, we specify a discretization, a macroscopic strain and a solver.

````@example custommaterial
disc = MoulinetSuquetDiscretization(microstructure)
macro_strain = MacroscopicStrain([0.01, 0, 0, 0, 0, 0])
solver = BasicScheme(mat_soft_custom.μ + mat_hard_custom.μ, microstructure)
sol = solve(microstructure, disc, macro_strain, solver)
nothing # hide
````

Now we repeat the same computation using the built-in model, which should give the exact same results.

````@example custommaterial
mat_soft_builtin = LinearIsotropicElastic{dim, Float64}(2100., 0.3)
mat_hard_builtin = LinearIsotropicElastic{dim, Float64}(72000., 0.22)
validation_microstructure = simple_microstructure(mat_soft_builtin, mat_hard_builtin)

disc = MoulinetSuquetDiscretization(validation_microstructure)
solver = BasicScheme(mat_soft_builtin.μ + mat_hard_builtin.μ, validation_microstructure)
validation_sol = solve(validation_microstructure, disc, macro_strain, solver)

sol.stress_avg ≈ validation_sol.stress_avg
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

