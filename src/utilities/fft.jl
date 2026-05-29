"""
    make_fft_plans(solver, field, field_fft, dim)

Create a forward `rfft` and inverse `irfft` plan for the spatial dimensions
`2:dim+1` of `field`, using the FFTW settings from `solver`.
"""
function make_fft_plans(
    solver::AbstractSolver,
    field::AbstractArray{T},
    field_fft::AbstractArray{Complex{T}},
    dim::Int
) where {T}
    FFTW.set_num_threads(solver.FFTW_num_threads)
    fft_plan  = plan_rfft(field, 2:dim+1; flags=solver.FFTW_flags)
    ifft_plan = plan_irfft(field_fft, size(field, 2), 2:dim+1; flags=solver.FFTW_flags)
    return fft_plan, ifft_plan
end