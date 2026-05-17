# load standard packages
using Pkg

# change into directory of this file
cd(@__DIR__)

# install non-registered package
Pkg.add(url="https://github.com/LudwigBoess/EntityTools.jl")
Pkg.add(url="https://github.com/LudwigBoess/MakieHelper.jl")

# activate environment and install packages
Pkg.activate(".")
Pkg.instantiate()