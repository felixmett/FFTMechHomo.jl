struct MicroStructure{dim, T <: AbstractFloat, A <: AbstractArray{<:AbstractMaterial}}
    materials::A       
end

function MicroStructure(materials::AbstractArray{<:AbstractMaterial})
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
            but materials are of spatial dimension $(first(dims))"
        )
    )
    MicroStructure{dim, T, typeof(materials)}(materials)
end

Base.ndims(::MicroStructure{dim, T, A}) where {dim, T <: AbstractFloat, A <: AbstractArray{<:AbstractMaterial}} = d
Base.size(microstructure::MicroStructure{dim, T, A}) where {dim, T <: AbstractFloat, A <: AbstractArray{<:AbstractMaterial}} = size(microstructure.materials)