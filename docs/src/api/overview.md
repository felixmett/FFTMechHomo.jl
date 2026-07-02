# Overview
`FFTMechHomo.jl` is organized around the following core components, which reflect the computational workflow of FFT-based homogenization:

- [Material](@ref): Constitutive laws for individual phases  
- [Microstructure](@ref): Spatial arrangement of material phases in the unit cell  
- [Discretization](@ref): Numerical representation of fields on the grid  
- [Solver](@ref): Algorithms for solving the Lippmann–Schwinger equation  
- [Loading](@ref): Definition of the applied macroscopic strain
