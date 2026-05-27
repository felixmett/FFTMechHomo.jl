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

function cell_average!(
    cell_avg::Vector{T}, 
    field_fft::AbstractArray{Complex{T}}, 
    consts::GridConstants{dim, T}
) where {dim, T}
    cell_avg .= real(field_fft[1:consts.n_voigt, consts.zero_idx]) ./ consts.n_cells
    return
end