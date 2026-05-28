"""
    MaterialGroup{dim, M}

Internal struct.
Groups indices by concrete material type to enable compile-time dispatch in compute_stress_field!.
"""
struct MaterialGroup{dim, M <: AbstractMaterial}
    indices::Vector{CartesianIndex{dim}}
end

"""
    InternalMicrostructure{dim, T, G, A}

Internal preprocessed representation of a [`Microstructure`](@ref).

Materials are grouped by concrete type into a `Tuple` of [`MaterialGroup`](@ref)s.
This allows `compute_stress_field!` to resolve `compute_stress!` dispatch once
per material type rather than once per voxel via `@generated` code emission.
"""
struct InternalMicrostructure{dim, T <: AbstractFloat, G <: Tuple, A <: AbstractArray{<: AbstractMaterial}}
    groups::G
    materials::A
end

"""
    InternalMicrostructure(ms::Microstructure)

Preprocess `ms` into an [`InternalMicrostructure`](@ref) by grouping voxel
indices by concrete material type.

This is called automatically by [`solve`](@ref) and does not need to be called
manually.
"""
function InternalMicrostructure(ms::Microstructure{dim, T, A}) where {dim, T <: AbstractFloat, A <: AbstractArray{<:AbstractMaterial}}
    type_map = Dict{Type, Vector{CartesianIndex{dim}}}()
    for i in CartesianIndices(ms.materials)
        mat_type = typeof(ms.materials[i])
        push!(get!(type_map, mat_type, Vector{CartesianIndex{dim}}()), i)
    end

    groups = Tuple(
        MaterialGroup{dim, mat_type}(indices)
        for (mat_type, indices) in type_map
    )

    InternalMicrostructure{dim, T, typeof(groups), typeof(ms.materials)}(groups, ms.materials)
end

Base.ndims(::InternalMicrostructure{dim}) where {dim} = dim
Base.size(microstructure::InternalMicrostructure) = size(microstructure.materials)

"""
    compute_stress_field!(stress, strain, ms::InternalMicrostructure)

Compute the stress response at every grid point of `ms` in-place.

# Extended help
This function is `@generated`. For a given [`InternalMicrostructure`](@ref)
type, it emits one explicit [`compute_group_stress!`](@ref) call per material group.
For an `InternalMicrostructure` with two groups this compiles to:

```julia
compute_group_stress!(stress, strain, ms.groups[1], ms.materials)
compute_group_stress!(stress, strain, ms.groups[2], ms.materials)
```

Each call carries a concrete `MaterialGroup{dim, M}` type, so `compute_stress!`
specialization happens at compile time. Subsequent calls reuse the compiled code.
"""
@generated function compute_stress_field!(
    stress::AbstractArray{T},
    strain::AbstractArray{T},
    ms::InternalMicrostructure{dim, T, G, A}
) where {dim, T, G, A}
    n = fieldcount(fieldtype(ms, :groups))
    calls = [:(compute_group_stress!(stress, strain, ms.groups[$i], ms.materials)) for i in 1:n]
    quote $(calls...) end
end

"""
    compute_group_stress!(stress, strain, group::MaterialGroup, materials::AbstractArray{<:AbstractMaterial})

Computes the stress response for each distinct material group inside of the materials of a microstructure.

!!! Note
    The ::M type assertion on materials[i] narrows the abstract array element type to the concrete `M` from
    `MaterialGroup{dim, M}`, allowing the compiler to resolve compute_stress! dispatch at compile time rather
    than run time.
"""
function compute_group_stress!(
    stress::AbstractArray{T},
    strain::AbstractArray{T},
    group::MaterialGroup{dim, M},
    materials::AbstractArray{<:AbstractMaterial}
) where {T, dim, M <: AbstractMaterial}
    for i in group.indices
        compute_stress!(stress, strain, materials[i]::M, i)
    end
end

"""
    compute_polarization!(stress, strain, ms::InternalMicrostructure, ref::ReferenceMaterial)

Compute the polarization field `τ = σ(ε) - σ₀(ε)` in-place, storing the result in `stress`.

On return, `stress` holds the polarization field rather than the physical stress.
Combines [`compute_stress_field!`](@ref) with [`subtract_from_stress_field!`](@ref) without requiring an additional buffer.
"""
function compute_polarization_field!(
    stress::AbstractArray{T},
    strain::AbstractArray{T},
    ms::InternalMicrostructure{dim, T},
    ref::ReferenceMaterial{dim, T}
) where {dim, T}
    compute_stress_field!(stress, strain, ms)
    subtract_from_stress_field!(stress, strain, ref)
    return
end