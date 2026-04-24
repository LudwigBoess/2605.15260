using DelimitedFiles
using CairoMakie, MakieHelper
using Colors

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

    shalaby = [1.0, 1.0, 1.18, 1.62, 1.83, 2.28]
    lines!(ax, 0:length(shalaby)-1, shalaby, color=shalabycolor, linewidth=2.5, label="Shalaby et al. 2017")

    lines!(ax, 0:length(x_labels)-1, deposit_timing, color=depositcolor, linewidth=2.5, label="Deposit")
    lines!(ax, 0:length(x_labels)-1, pusher_timing, color=pushercolor, linewidth=2.5, label="Pusher")
    lines!(ax, 0:length(x_labels)-1, total_timing, color=totalcolor, linewidth=4.5, label="Total")

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

set_light_theme!(fontsize=20)

data_path = data_path * "performance/"
plot_name = plot_path * "Fig20.pdf"

Fig20(data_path, plot_name)