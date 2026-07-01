#! /bin/bash
julia --project=docs/ -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.update()'

julia --project=docs/ docs/literatedocs.jl

julia --project=docs/ docs/make.jl