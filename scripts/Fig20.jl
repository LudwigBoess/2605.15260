#include(joinpath(@__DIR__, "config.jl"))

using DelimitedFiles
using Printf
using CairoMakie, MakieHelper
using Colors
using EntityTools

"""
    parse_timing(run_dir, output_file; skip_steps=5)

Walk through every run subfolder in `run_dir`, parse its `weibel/weibel.out`
with `EntityTools.parse_output`, and write the per-method average timings
(total, Communications, CurrentDeposit, ParticlePusher) in ms to
`output_file`, skipping the first `skip_steps` timesteps.
"""
function parse_timing(run_dir::AbstractString, output_file::AbstractString;
                      skip_steps::Integer = 5)

    method_dirs = ["zigzag",
                     "esirkepov_1st", "esirkepov_2nd", "esirkepov_3rd",
                     "esirkepov_4th", "esirkepov_5th", "esirkepov_6th",
                     "esirkepov_7th", "esirkepov_8th", "esirkepov_9th",
                     "esirkepov_10th", "esirkepov_11th"]

    method_labels = ["zigzag",
                       "O1", "O2", "O3", "O4", "O5", "O6",
                       "O7", "O8", "O9", "O10", "O11"]

    open(output_file, "w") do io
        println(io, "method; total; comm; deposit; pusher")
        for (sub, label) in zip(method_dirs, method_labels)
            outpath = joinpath(run_dir, sub, "weibel", "weibel.out")
            if !isfile(outpath)
                @warn "missing output file" outpath
                continue
            end

            sim   = parse_output(outpath)
            steps = sim.steps[skip_steps+1:end]
            n     = length(steps)

            tot  = 1e3 * sum(s.total                       for s in steps) / n
            comm = 1e3 * sum(s.substeps["Communications"]  for s in steps) / n
            dep  = 1e3 * sum(s.substeps["CurrentDeposit"]  for s in steps) / n
            pus  = 1e3 * sum(s.substeps["ParticlePusher"]  for s in steps) / n

            @printf(io, "%s; %.5e; %.5e; %.5e; %.5e\n", label, tot, comm, dep, pus)
        end
    end
    return output_file
end

function read_timing(fi)

    d = readdlm(fi, ';', header=true)[1]

    deposit_timing = d[:,4]
    deposit_timing ./= deposit_timing[1]
    pusher_timing = d[:,5]
    pusher_timing ./= pusher_timing[1]
    total_timing = d[:,2]
    total_timing ./= total_timing[1]

    return deposit_timing, pusher_timing, total_timing
end

function Fig20(data_path, plot_name)

    x_labels = [L"\mathcal{Z}", L"\mathcal{S}_1", 
                L"\mathcal{S}_2", L"\mathcal{S}_3", 
                L"\mathcal{S}_4", L"\mathcal{S}_5", 
                L"\mathcal{S}_6", L"\mathcal{S}_7", 
                L"\mathcal{S}_8", L"\mathcal{S}_9",
                L"\mathcal{S}_{10}", L"\mathcal{S}_{11}"]

    depositcolor = RGBA(0.1725, 0.6275, 0.1725, 1.0)
    pushercolor = RGBA(1.0, 0.4980, 0.0549, 1.0)
    totalcolor = RGBA(0.1216, 0.4667, 0.7059, 1.0)
    shalabycolor = RGBA(0.8392, 0.1529, 0.1569, 1.0)

    f = Figure(size=(1300,400))
    ax = Axis(f[1,1], yscale=log10, 
                xlabel="Interpolation method", ylabel="Runtime increase",
              yminorticks=IntervalsBetween(9),
              xticks=0:length(x_labels)-1,
              xtickformat = x -> x_labels[Int64.(x).+1],
              xminorticksvisible=false,
              title = "1D"
              )
    ylims!(0.1, 1.0e3)
    xlims!(0.0, length(x_labels)-1)
    ax.xminorticks = false

    deposit_timing, pusher_timing, total_timing = read_timing(data_path * "1D.txt")


    lines!(ax, 0:length(x_labels)-1, deposit_timing, color=depositcolor, linewidth=2.5, label="Deposit")
    lines!(ax, 0:length(x_labels)-1, pusher_timing, color=pushercolor, linewidth=2.5, label="Pusher")
    lines!(ax, 0:length(x_labels)-1, total_timing, color=totalcolor, linewidth=4.5, label="Total")
    
    # renormalize to our 1st order for direct comparison
    shalaby = [1.0, 1.18, 1.62, 1.83, 2.28] .* total_timing[2]
    lines!(ax, 1:length(shalaby), shalaby, color=shalabycolor, linewidth=2.5, label="Shalaby et al. 2017")

    axislegend(ax, framevisible=false, position = :lt)

    ax = Axis(f[1,2], yscale=log10, 
                xlabel="Interpolation method", ylabel="Runtime increase",
              yminorticks=IntervalsBetween(9),
              xticks=0:length(x_labels)-1,
              xtickformat = x -> x_labels[Int64.(x).+1],
              xminorticksvisible=false,
              title = "2D"
              )
    ylims!(0.1, 1.0e3)
    xlims!(0.0, length(x_labels)-1)
    ax.xminorticks = false

    deposit_timing, pusher_timing, total_timing = read_timing(data_path * "2D.txt")

    lines!(ax, 0:length(x_labels)-1, deposit_timing, color=depositcolor, linewidth=2.5, label="Deposit")
    lines!(ax, 0:length(x_labels)-1, pusher_timing, color=pushercolor, linewidth=2.5, label="Pusher")
    lines!(ax, 0:length(x_labels)-1, total_timing, color=totalcolor, linewidth=4.5, label="Total")

    ax = Axis(f[1,3], yscale=log10, 
                xlabel="Interpolation method", ylabel="Runtime increase",
              yminorticks=IntervalsBetween(9),
              xticks=0:length(x_labels)-1,
              xtickformat = x -> x_labels[Int64.(x).+1],
              xminorticksvisible=false,
              title = "3D"
              )
    ylims!(0.1, 1.0e3)
    xlims!(0.0, length(x_labels)-1)

    deposit_timing, pusher_timing, total_timing = read_timing(data_path * "3D.txt")

    lines!(ax, 0:length(x_labels)-1, deposit_timing, color=depositcolor, linewidth=2.5, label="Deposit")
    lines!(ax, 0:length(x_labels)-1, pusher_timing, color=pushercolor, linewidth=2.5, label="Pusher")
    lines!(ax, 0:length(x_labels)-1, total_timing, color=totalcolor, linewidth=4.5, label="Total")

    save(plot_name, f)
end

set_dark_theme!(fontsize=20)

perf_in_path  = "/lus/flare/projects/CRsInCosmoSims/PIC/Paper/entity_convergence/performance"
perf_out_path = "/lus/flare/projects/CRsInCosmoSims/PIC/Paper/2600.0000/simulations/performance"
plot_path = "/lus/flare/projects/CRsInCosmoSims/PIC/Paper/2600.0000/Plots/"

#for dim in ("1D", "2D", "3D")
#    parse_timing(joinpath(perf_in_path,  dim),
#                 joinpath(perf_out_path, dim * ".txt"), skip_steps=2)
#end

data_path = perf_out_path * "/"
plot_name = plot_path * "Fig20.png"

Fig20(data_path, plot_name)