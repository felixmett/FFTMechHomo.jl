struct BasicScheme{T <: AbstractFloat} <: AbstractSolver
    α₀::T
    tol::T
    maxiter::Int
    FFTW_flags::UInt32
    FFTW_num_threads::Int
end

BasicScheme(
    α₀::T; tol=T(1e-5), maxiter=100,
    FFTW_flags=FFTW.MEASURE, FFTW_num_threads=1
) where {T <: AbstractFloat} = BasicScheme{T}(α₀, tol, maxiter, FFTW_flags, FFTW_num_threads)