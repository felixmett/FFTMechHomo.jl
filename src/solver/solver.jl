abstract type AbstractSolver end
abstract type LinearSolver <: AbstractSolver end
abstract type NonLinearSolver <: AbstractSolver end

# fft_plan  = plan_rfft(τ, 2:d+1; flags=FFTW.MEASURE)
# ifft_plan = plan_irfft(τ_hat, size(ε, 2), 2:d+1; flags=FFTW.MEASURE)
# why not fft_plan! ???

# # Somehow make this aware if linear or not
# struct Solution{S <: AbstractSolver}
#     stress
#     strain
#     stress_avg
#     residuals
# end

# function residual!(
#     σ::AbstractArray{T},
#     σ₀::AbstractArray{T},   # reuse as scratch
#     ε::AbstractArray{T},
#     ε_prev::AbstractArray{T},
#     σ_bar::AbstractVector{T},
#     ref::ReferenceMaterial
# ) where {T}
#     # σ₀ reused as scratch: compute stress of (ε - ε_prev) inplace
#     @. σ₀ = ε - ε_prev
#     compute_stress!(σ, ref, σ₀)
#     norm_factor = T(prod(size(ε)[2:end]))
#     nom   = sqrt((1/norm_factor) * dot(σ, σ))
#     denom = sqrt(dot(σ_bar, σ_bar))
#     return nom / denom
# end