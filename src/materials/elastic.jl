abstract type AbstractElastic <: AbstractMaterial end

struct LinearIsotropicElastic{T <: AbstractFloat} <: AbstractElastic
    E::T
    nu::T

    function LinearIsotropicElastic(E::T, nu::T) where T <: AbstractFloat
        E > 0 || throw(ArgumentError("E must be positive"))
        -1 < nu < 0.5 || throw(ArgumentError("nu must be in (-1, 0.5)"))
        new{T}(E, nu)
    end
end