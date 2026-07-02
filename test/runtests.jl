using FFTMechHomo
using FFTW
using LinearAlgebra
using Test

include("discretization/greenoperator.jl")
include("discretization/moulinetsuquet.jl")

include("material/models/elastic.jl")
include("material/models/reference.jl")

include("microstructure/internal.jl")
include("microstructure/microstructure.jl")

include("solver/basic.jl")
include("solver/macroscopicstrain.jl")
include("solver/solver.jl")

include("utilities/tensors.jl")

include("integration.jl")