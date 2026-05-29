struct BasicScheme{L <: AbstractLinearity, T <: AbstractFloat} <: AbstractSolver
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
    L = Linear # TODO: make that inferrable from a microstructure!
    return BasicScheme{L, T}(T(α₀), T(tol), maxiter, FFTW_flags, FFTW_num_threads)
end

function solve(
    microstructure::Microstructure{dim, T},
    disc::AbstractDiscreteGreenOperator{dim, T},
    macro_strain::MacroscopicStrain{dim, T},
    solver::BasicScheme{Linear, T}
) where {dim, T}
    # reference material and internal microstructure
    ref = ReferenceMaterial{dim, T}(solver.α₀)
    internal_microstructure = InternalMicrostructure(microstructure)
    grid_constants = GridConstants(internal_microstructure)

    # preallocate all buffers
    strain = initialize_strain_field(macro_strain, internal_microstructure)  # (n_voigt, nx, ...)
    strain_prev = similar(strain)
    polarization = similar(strain)
    polarization_fft   = Array{Complex{T}}(undef, rfft_output_size(disc))

    # stress from macroscopic strain, used in residual denominator
    macro_stress_ref = compute_stress(macro_strain.data, ref)
    stress_avg = similar(macro_stress_ref)

    fft_plan, ifft_plan = make_fft_plans(solver, strain, polarization_fft, dim)

    res = Vector{T}(undef, solver.maxiter + 1)
    res[1] = T(Inf)
    k = 0

    while k < solver.maxiter && res[k+1] > solver.tol
        k += 1

        # save previous strain for residual
        copyto!(strain_prev, strain)

        # compute polarization τ = σ(ε) - σ₀(ε) inplace, no alloc
        compute_polarization_field!(polarization, strain, internal_microstructure, ref)

        # forward FFT of polarization
        mul!(polarization_fft, fft_plan, polarization)

        # average stress: σ_bar = σ_δε + Re(τ_hat[:,1,1,...]) / norm_factor
        cell_average!(stress_avg, polarization_fft, grid_constants)
        stress_avg .+= macro_stress_ref

        # apply Green operator inplace
        Γ⁰!(polarization_fft, ref, disc)

        # impose macroscopic strain at zero frequency
        impose_macroscopic_strain!(polarization_fft, macro_strain, grid_constants)

        # inverse FFT back into ε
        mul!(strain, ifft_plan, polarization_fft)

        res[k+1] = residual!(polarization, strain, strain_prev, stress_avg, ref, grid_constants)
    end

    # final stress
    compute_stress_field!(polarization, strain, internal_microstructure)
    mean!(stress_avg, polarization)

    return Solution{T}(strain, stress, stress_avg, res[2:k+1], k < solver.maxiter, k)
end