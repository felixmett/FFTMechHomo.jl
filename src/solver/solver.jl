"""
    AbstractSolver

Abstract supertype for all FFT-based homogenization solvers.

# Interface
Every solver must implement:
- `solve(solver, microstructure, ref, disc, bc)`
"""
abstract type AbstractSolver end

"""
    solve(microstructure, discretization, macro_strain, solver)

Solve the periodic homogenization problem for a given macroscopic strain.

Returns a [`Solution`](@ref) containing the strain and stress fields, the
macroscopic stress average, convergence history, and iteration count.

Arguments:
- `microstructure`: spatial distribution of material models
- `discretization`: discretization of the Green operator in Fourier space
- `macro_strain`: applied macroscopic strain
- `solver`: controls the solution scheme, convergence tolerance, and FFT settings
"""
function solve end

"""
    AbstractLinearity

Abstract supertype for linearity markers used as type parameters in solver structs.
See [`Linear`](@ref) and [`NonLinear`](@ref).
"""
abstract type AbstractLinearity end

"""
    Linear <: AbstractLinearity

Marks a solver as operating on a microstructure with only [`HistoryIndependent`](@ref)
materials, enabling a single-loop solution scheme.
"""
abstract type Linear <: AbstractLinearity end

"""
    NonLinear <: AbstractLinearity

Marks a solver as operating on a microstructure containing [`HistoryDependent`](@ref)
materials, requiring an inner iteration loop.
"""
abstract type NonLinear <: AbstractLinearity end

struct GridConstants{dim, T <: AbstractFloat}
    n_voigt::Int
    n_cells::T
    zero_idx::CartesianIndex{dim}
end

"""
    GridConstants(microstructure)

Internally used structure that precomputes and holds grid constants used throughout the solver loop.
"""
function GridConstants(microstructure::InternalMicrostructure{dim, T}) where {dim, T}
   n_voigt = dim^2 - dim^(dim-2)
   n_cells = T(prod(size(microstructure.materials)))
   zero_idx = CartesianIndex(ntuple(x -> 1, dim))
   return GridConstants{dim, T}(n_voigt, n_cells, zero_idx)
end

"""
    residual!(stress, strain, strain_prev, stress_avg, ref, consts)

Compute the normalized stress residual between the current and previous strain iterate.

Overwrites `stress` with the reference stress of `strain - strain_prev`.

!!! warning
    `stress` is overwritten inside of this function.
"""
function residual!(
    stress::AbstractArray{T},
    strain::AbstractArray{T},
    strain_prev::AbstractArray{T},
    stress_avg::AbstractVector{T},
    ref::ReferenceMaterial{dim, T},
    consts::GridConstants{dim, T}
) where {dim, T}
    stress .= strain .- strain_prev
    compute_stress_field!(stress, stress, ref)

    inner_product_stress = inner_product(stress, consts.n_cells, :stress)
    norm_stress_avg = tensor_dot(stress_avg, :stress)
    return sqrt(inner_product_stress / norm_stress_avg)
end

"""
    cell_average!(cell_avg, field_fft, consts)

Compute the cell average of a real-valued field from its Fourier representation.
Overwrites `cell_avg` to store the result.
"""
function cell_average!(
    cell_avg::Vector{T},
    field_fft::AbstractArray{Complex{T}},
    consts::GridConstants{dim, T}
) where {dim, T}
    cell_avg .= real(field_fft[1:consts.n_voigt, consts.zero_idx]) ./ consts.n_cells
    return
end

"""
    impose_macroscopic_strain!(field_fft, macro_strain, consts)

Enforce the macroscopic strain constraint `⟨ε⟩ = ε̄` by setting the zero-frequency mode of `field_fft` in-place.

In the periodic homogenization problem the strain field is decomposed into the applied macroscopic strain and
a zero-mean periodic fluctuation. The zero-frequency mode of the DFT corresponds to the spatial average,
so setting it to the macroscopic strain enforces the mean strain constraint after the inverse FFT.
"""
function impose_macroscopic_strain!(
    field_fft::AbstractArray{Complex{T}},
    macro_strain::MacroscopicStrain{dim, T},
    consts::GridConstants{dim, T}
) where {dim, T}
    field_fft[1:consts.n_voigt, consts.zero_idx] .= Complex{T}.(macro_strain.data .* consts.n_cells)
    return
end