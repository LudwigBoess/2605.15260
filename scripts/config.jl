# path to the package envirnoment
const global dir_path = joinpath(@__DIR__, "..")

# activate the environment
using Pkg;
Pkg.activate(dir_path);


# paths to data, plots and maps
const global data_path = dir_path * "/simulations/"
const global plot_path = dir_path * "/Plots/"
