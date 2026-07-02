"""
    BasicScheme(α₀, microstructure; tol, maxiter, FFTW_flags, FFTW_num_threads)

Moulinec-Suquet basic scheme solver for periodic homogenization.

The reference material parameter `α₀` controls the convergence of the solution scheme.

# Arguments
- `α₀`: reference material parameter, see [`ReferenceMaterial`](@ref)
- `microstructure`: used to infer the numeric type `T` and the solver's linearity marker `L`

# Keyword Arguments
- `tol`: convergence tolerance on the stress residual (default: `1e-5`)
- `maxiter`: maximum number of iterations (default: `100`)
- `FFTW_flags`: FFTW planner flags (default: `FFTW.MEASURE`)
- `FFTW_num_threads`: number of threads for FFTW (default: `1`)

# Example
```julia
ms = Microstructure(fill(LinearIsotropicElastic{3}(72e3, 0.3), 32, 32, 32))
solver = BasicScheme(2 * 72e3, ms; tol=1e-6, maxiter=200)
```
"""
struct BasicScheme{T <: AbstractFloat} <: AbstractSolver
    α₀::T
    tol::T
    maxiter::Int
    FFTW_flags::UInt32
    FFTW_num_threads::Int
end

function BasicScheme(
    α₀::Real, microstructure::Microstructure{dim, T};
    tol=1e-5, maxiter=100, FFTW_flags=FFTW.MEASURE, FFTW_num_threads=1
) where {dim, T <: AbstractFloat}
    return BasicScheme{T}(T(α₀), T(tol), maxiter, FFTW_flags, FFTW_num_threads)
end

"""
    solve(microstructure, disc, macro_strain, solver::BasicScheme)

Solve the periodic homogenization problem using the Moulinec-Suquet basic scheme.

See [`solve`](@ref) for the general interface and [`Solution`](@ref) for the
return type.
"""
function solve(
    microstructure::Microstructure{dim, T},
    disc::AbstractDiscreteGreenOperator{dim, T},
    macro_strain::MacroscopicStrain{dim, T},
    solver::BasicScheme{T}
) where {dim, T}
    ref = ReferenceMaterial{dim, T}(solver.α₀)
    internal_microstructure = InternalMicrostructure(microstructure)
    grid_constants = GridConstants(internal_microstructure)
    macro_stress_ref = compute_stress(macro_strain.data, ref)

    # preallocate all buffers
    strain = initialize_strain_field(macro_strain, internal_microstructure)  # (n_voigt, nx, ...)
    strain_prev = similar(strain)
    polarization = similar(strain)
    polarization_fft   = Array{Complex{T}}(undef, rfft_output_size(disc))
    stress_avg = similar(macro_stress_ref)

    fft_plan, ifft_plan = make_fft_plans(solver, strain, polarization_fft, dim)
    # make sure to not use intermediate garbage from FFT plan creation
    strain .= initialize_strain_field(macro_strain, internal_microstructure)

    res = Vector{T}(undef, solver.maxiter + 1)
    res[1] = T(Inf)
    k = 0
    # Buffer roles throughout the iteration loop:
    #   strain:             current strain field, updated via IFFT of polarization_fft each iteration
    #   strain_prev:        strain from previous iteration, used in residual
    #   polarization:       τ = σ(ε) - σ₀(ε),
    #                       note that this is also used as a buffer for stresses σ
    #                       inside of the residual computation and for the final stress computation
    #                       in the end
    #   polarization_fft:   Fourier transform of the polarization τ,
    #                       note that after the in-place green operator `Γ⁰!` is applied, the field
    #                       holds the current value of the Fourier transformed strain
    #   stress_avg:         macroscopic stress ⟨σ⟩, recomputed each iteration
    while k < solver.maxiter && res[k+1] > solver.tol
        k += 1
        copyto!(strain_prev, strain)

        compute_polarization_field!(polarization, strain, internal_microstructure, ref)
        mul!(polarization_fft, fft_plan, polarization)

        cell_average!(stress_avg, polarization_fft, grid_constants)
        stress_avg .+= macro_stress_ref

        Γ⁰!(polarization_fft, ref, disc)
        impose_macroscopic_strain!(polarization_fft, macro_strain, grid_constants)
        mul!(strain, ifft_plan, polarization_fft)

        res[k+1] = residual!(polarization, strain, strain_prev, stress_avg, ref, grid_constants)
    end

    # final stress
    compute_stress_field!(polarization, strain, internal_microstructure)
    mean!(stress_avg, polarization)

    return Solution{T}(strain, polarization, stress_avg, res[2:k+1], k < solver.maxiter, k)
end