module FFTMechHomo

using FFTW
using LinearAlgebra
using StaticArrays
using Statistics

# types
export HistoryIndependent

# structs
export MoulinetSuquetDiscretization
export LinearIsotropicElastic
export Microstructure
export BasicScheme
export MacroscopicStrain

# functions
export solve
export compute_stress!

include("material/material.jl")
include("microstructure/microstructure.jl")
include("discretization/greenoperator.jl")
include("material/models/elastic.jl")
include("material/models/reference.jl")
include("microstructure/internal.jl")
include("solver/macroscopicstrain.jl")
include("utilities/tensors.jl")
include("solver/solution.jl")
include("solver/solver.jl")
include("utilities/fft.jl")

include("discretization/moulinetsuquet.jl")

include("solver/basic.jl")



end


