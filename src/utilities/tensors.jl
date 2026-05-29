"""
    voigt_weights(::Val{N}, ::Type{T}, ::Val{kind})

Return the diagonal weight matrix for Voigt notation tensor dot products in 2D and 3D.

For a symmetric second-order tensor stored in Voigt notation, the correct inner product requires weights
to account for the off-diagonal components appearing only once in the vector but twice in the full tensor
"""
voigt_weights(::Val{3}, ::Type{T}, ::Val{:strain}) where {T} = Diagonal(T[1, 1, 0.5])
voigt_weights(::Val{6}, ::Type{T}, ::Val{:strain}) where {T} = Diagonal(T[1, 1, 1, 0.5, 0.5, 0.5])
voigt_weights(::Val{3}, ::Type{T}, ::Val{:stress}) where {T} = Diagonal(T[1, 1, 2])
voigt_weights(::Val{6}, ::Type{T}, ::Val{:stress}) where {T} = Diagonal(T[1, 1, 1, 2, 2, 2])
voigt_weights(::Val{3}, ::Type{T}, ::Val{:mixed}) where {T} = Diagonal(T[1, 1, 1])
voigt_weights(::Val{6}, ::Type{T}, ::Val{:mixed}) where {T} = Diagonal(T[1, 1, 1, 1, 1, 1])
voigt_weights(::Val{N}, ::Type{T}, ::Val{S}) where {N, T, S <: Symbol} = throw(ArgumentError("Vectors must be of length 3 or 6, got $N"))

"""
    tensor_dot(u, v, kind)
    tensor_dot(u, kind)

Compute the weighted inner product of two Voigt vectors `u` and `v`, such that it represents
the tensor dot products of their tensor equivalents `U` and `V`:

    U : V = uᵀ * W * v

where `W` is the Voigt weight matrix for `kind`:
- `:strain` — weights for strain-strain inner product
- `:stress` — weights for stress-stress inner product
- `:mixed`  — identity weights for stress-strain inner product (two-argument form only)

See [`voigt_weights`](@ref) for the explicit weight values.

The single-argument form computes `U : U` and is only meaningful for `:strain` and `:stress`.
"""
function tensor_dot(u::AbstractVector{T}, v::AbstractVector{T}, kind::Symbol) where {T <: AbstractFloat}
    length(u) == length(v) || throw(ArgumentError("Vectors must have equal length, got $(length(u)) and $(length(v))"))
    W = voigt_weights(Val(length(u)), T, Val(kind))
    return transpose(u) * W * v
end

tensor_dot(u::AbstractVector{T}, kind::Symbol) where {T <: AbstractFloat} = tensor_dot(u, u, kind)

"""
    inner_product(U, V, volume, kind)

Compute the volume-averaged inner product of two tensor fields:

    (U, V) = 1/|Q| ∫_Q U : V dx

where `:` denotes the tensor dot product computed in Voigt notation through appropriate weights determined by `kind`.
See [`tensor_dot`](@ref) for the pointwise operation.
"""
function inner_product(
    U::AbstractArray{T},
    V::AbstractArray{T},
    volume::T,
    kind::Symbol
) where {T}
    out = T(0)
    for i in CartesianIndices(size(U)[2:end])
        out += tensor_dot(view(U, :, i), view(V, :, i), kind)
    end
    return out / volume
end

inner_product(U::AbstractArray{T}, volume::T, kind::Symbol) where {T} = inner_product(U, U, volume, kind)