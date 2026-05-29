"""
    Solution

Result returned by [`solve`](@ref).

# Fields
- `strain`: converged strain field
- `stress`: converged stress field
- `stress_avg`: macroscopic stress average `⟨σ⟩`
- `residuals`: residual norms at each iteration
- `converged`: `true` if the solver converged within `maxiter` iterations
- `iterations`: number of iterations performed
"""
struct Solution{T}
    strain::Array{T}
    stress::Array{T}
    stress_avg::Vector{T}
    residuals::Vector{T}
    converged::Bool
    iterations::Int
end