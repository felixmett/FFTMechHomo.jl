voigt_weights(::Val{3}, ::Type{T}, ::Val{:strain}) where {T} = Diagonal(T[1, 1, 0.5])
voigt_weights(::Val{6}, ::Type{T}, ::Val{:strain}) where {T} = Diagonal(T[1, 1, 1, 0.5, 0.5, 0.5])
voigt_weights(::Val{3}, ::Type{T}, ::Val{:stress}) where {T} = Diagonal(T[1, 1, 2])
voigt_weights(::Val{6}, ::Type{T}, ::Val{:stress}) where {T} = Diagonal(T[1, 1, 1, 2, 2, 2])
voigt_weights(::Val{N}, ::Type{T}, ::Val{S}) where {N, T, S <: Symbol} = throw(ArgumentError("Vectors must be of length 3 or 6, got $N"))

function tensor_dot(u::AbstractVector{T}, v::AbstractVector{T}, kind::Symbol) where {T <: AbstractFloat}
    length(u) == length(v) || throw(ArgumentError("Vectors must have equal length, got $(length(u)) and $(length(v))"))
    W = voigt_weights(Val(length(u)), T, Val(kind))
    return transpose(u) * W * v
end

tensor_dot(u::AbstractVector{T}, kind::Symbol) where {T <: AbstractFloat} = tensor_dot(u, u, kind)