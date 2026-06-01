struct Microstructure{dim, T <: AbstractFloat, A <: AbstractArray{<:AbstractMaterial}}
    materials::A
end

"""
    Microstructure(materials)

Construct a microstructure from a spatial array of material models.

Each element of `materials` corresponds to one voxel of the microstructure.
All materials must have consistent spatial dimension and numeric type.

# Arguments
- `materials`: array of [`AbstractMaterial`](@ref) instances with shape `(n₁, n₂)` for 2D or `(n₁, n₂, n₃)` for 3D

# Example
```julia
mat = LinearIsotropicElastic{3}(210e3, 0.3)
ms  = Microstructure(fill(mat, 32, 32, 32))
```
"""
function Microstructure(materials::AbstractArray{<:AbstractMaterial})
    dim = ndims(first(materials))
    T = eltype(first(materials))

    for mat in materials
        ndims(mat) == dim || throw(ArgumentError(
            "inconsistent material dimensions: expected $dim, got $(ndims(mat))"
        ))
        eltype(mat) == T || throw(ArgumentError(
            "inconsistent material types: expected $T, got $(eltype(mat))"
        ))
    end

    length(size(materials)) == dim || throw(
        ArgumentError(
            "microstructure has $(ndims(materials)) spatial dimensions,
            but materials are of spatial dimension $dim"
        )
    )
    Microstructure{dim, T, typeof(materials)}(materials)
end

Base.ndims(::Microstructure{dim, T, A}) where {dim, T <: AbstractFloat, A <: AbstractArray{<:AbstractMaterial}} = dim
Base.size(microstructure::Microstructure{dim, T, A}) where {dim, T <: AbstractFloat, A <: AbstractArray{<:AbstractMaterial}} = size(microstructure.materials)