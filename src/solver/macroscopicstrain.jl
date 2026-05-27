abstract type AbstractPrescribedStrain end

struct MacroscopicStrain{dim, T <: AbstractFloat} <: AbstractPrescribedStrain
    data::Vector{T}
end

struct MacroscopicStrainSteps{dim, T <: AbstractFloat} <: AbstractPrescribedStrain
    steps::Vector{Vector{T}}  # one strain vector per load step
end

initialize_strain_field(
    bc::MacroscopicStrain{dim, T},
    ms::InternalMicrostructure{dim, T}
) where {dim, T} = repeat(bc.data, 1, size(ms.materials)...)

function impose_macroscopic_strain!(
    field_fft::AbstractArray{Complex{T}},
    bc::MacroscopicStrain{dim, T},
    consts::GridConstants{dim, T}
) where {dim, T}
    field_fft[1:consts.n_voigt, consts.zero_idx] .= Complex{T}.(bc.data .* consts.n_cells)
    return
end 