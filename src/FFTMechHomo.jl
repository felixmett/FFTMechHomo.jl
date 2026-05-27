module FFTMechHomo

using FFTW
using LinearAlgebra
using StaticArrays
using Statistics

# structs
export MoulinetSuquetDiscretization
export LinearIsotropicElastic
export Microstructure
export BasicScheme
export MacroscopicStrain

# functions
export solve

include("material/material.jl")
include("microstructure/microstructure.jl")
include("discretization/greenoperator.jl")
include("material/models/elastic.jl")
include("material/models/reference.jl")
include("microstructure/internal.jl")
include("solver/solver.jl")

include("discretization/moulinetsuquet.jl")

include("solver/macroscopicstrain.jl")
include("solver/basic.jl")

include("utilities/fft.jl")
include("utilities/tensors.jl")

end


