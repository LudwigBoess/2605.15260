include(joinpath(@__DIR__, "config.jl"))

using CairoMakie, MakieHelper
using Printf
using DelimitedFiles
using EntityTools
using Statistics
using LaTeXStrings

function divE_err(sim_path)

    info = EntityData(sim_path)

    d = read_field(info, ["DivE"], t=[0.0, Inf], verbose=true, reduction_function=mean)

    return d["t"], sqrt.(d["DivE"].^2)
end

function write_divE_error(sim_paths, output_file)

    mean_divE = Vector{Float64}(undef, length(sim_paths))
    std_divE = Vector{Float64}(undef, length(sim_paths))
    min_divE = Vector{Float64}(undef, length(sim_paths))
    max_divE = Vector{Float64}(undef, length(sim_paths))

    for i = 1:length(sim_paths)

        println("Processing folder: ", sim_paths[i])

        t, divE = divE_err(sim_paths[i])

        divE = log10.(divE[abs.(divE) .> 0.0])

        println("mean: ", mean(divE), " std: ", std(divE), " min: ", minimum(divE), " max: ", maximum(divE))
        mean_divE[i] = mean(divE)
        std_divE[i] = std(divE)
        min_divE[i] = minimum(divE)
        max_divE[i] = maximum(divE)

    end

    data = hcat(mean_divE, std_divE, min_divE, max_divE)
    writedlm(output_file, data)

end


function Fig03(data_paths, sims, plot_name)

    f = Figure(size=(625, 650))

    ypos_label = 6.e-19
    lw_dash = 4.0

    ax = Axis(f[1, 1], aspect = 3,
                xticks=1:length(sims),
                xtickformat = x -> sims[Int64.(x)],
                yscale = log10,
                yticks = [1.e-20, 1.e-19, 1.e-18],
                yminorticks = IntervalsBetween(9),
                xminorticksvisible=false,
                ytickformat = values -> [iszero(value) ? "0" : rich("$(value < 0 ? "-" : "")10", superscript("$(@sprintf("%i", log10(abs(value))))")) for value in values]
                )
    xlims!(ax, (0.5, length(sims) + 0.5))
    ylims!(ax, (1.e-20, 1.e-18))
    hidexdecorations!(ax, ticks=false)

    data = readdlm(data_paths[1])
    mean_divE = data[:, 1]
    std_divE = data[:, 2]
    min_divE = data[:, 3]
    max_divE = data[:, 4]
    rho_mean = 10 .^ mean_divE
    rho_upper = rho_mean .- 10 .^ (mean_divE .- std_divE)
    rho_lower = 10 .^ (mean_divE .+ std_divE) .- rho_mean
    rho_max = 10 .^ max_divE

    lines!(ax, 1:length(sims), rho_max, color = 1:length(sims), colormap=:batlow, colorrange=(0,length(sims)+1), linestyle=:dash,
        linewidth=lw_dash)
    band!(ax, 1:length(sims), rho_lower, rho_upper, color=1:length(sims), colormap=:batlow, colorrange=(0,length(sims)+1), alpha=0.8)
    scatter!(ax, 1:length(sims), rho_mean, markersize = 15, strokecolor=:white, strokewidth = 1.0,
            color=1:length(sims), colormap = :batlow, colorrange=(0,length(sims)+1))

    textlabel!(ax, Point2f(length(sims)-0.5, ypos_label), text = "1D",
                    offset = (0, 0), text_color = :black, alpha=0.0, strokewidth=0)

    ax = Axis(f[2, 1], aspect = 3,
                ylabel=L"\overline{\nabla \cdot \mathbf{E}} \quad [4\pi (\rho_+ - \rho_-)]",
                xticks=1:length(sims),
                xtickformat = x -> sims[Int64.(x)],
                yscale = log10,
                yticks = [1.e-20, 1.e-19, 1.e-18],
                yminorticks = IntervalsBetween(9),
                xminorticksvisible=false,
                ytickformat = values -> [iszero(value) ? "0" : rich("$(value < 0 ? "-" : "")10", superscript("$(@sprintf("%i", log10(abs(value))))")) for value in values]
                )
    xlims!(ax, (0.5, length(sims) + 0.5))
    ylims!(ax, (1.e-20, 1.e-18))
    hidexdecorations!(ax, ticks=false)

    data = readdlm(data_paths[2])
    mean_divE = data[:, 1]
    std_divE = data[:, 2]
    min_divE = data[:, 3]
    max_divE = data[:, 4]
    rho_mean = 10 .^ mean_divE
    rho_upper = rho_mean .- 10 .^ (mean_divE .- std_divE)
    rho_lower = 10 .^ (mean_divE .+ std_divE) .- rho_mean

    lines!(ax, 1:length(sims), rho_max, color = 1:length(sims), colormap=:batlow, colorrange=(0,length(sims)+1), linestyle=:dash,
        linewidth=lw_dash)
    band!(ax, 1:length(sims), rho_lower, rho_upper, color=1:length(sims), colormap=:batlow, colorrange=(0,length(sims)+1), alpha=0.8)
    scatter!(ax, 1:length(sims), rho_mean, markersize = 15, strokecolor=:white, strokewidth = 1.0,
            color=1:length(sims), colormap = :batlow, colorrange=(0,length(sims)+1))
    textlabel!(ax, Point2f(length(sims)-0.5, ypos_label), text = "2D",
                    offset = (0, 0), text_color = :black, alpha=0.0, strokewidth=0)
    ax = Axis(f[3, 1], aspect = 3,
                xticks=1:length(sims),
                xtickformat = x -> sims[Int64.(x)],
                yscale = log10,
                yticks = [1.e-20, 1.e-19, 1.e-18],
                yminorticks = IntervalsBetween(9),
                xminorticksvisible=false,
                ytickformat = values -> [iszero(value) ? "0" : rich("$(value < 0 ? "-" : "")10", superscript("$(@sprintf("%i", log10(abs(value))))")) for value in values]
                )
    xlims!(ax, (0.5, length(sims) + 0.5))
    ylims!(ax, (1.e-20, 1.e-18))

    data = readdlm(data_paths[3])
    mean_divE = data[:, 1]
    std_divE = data[:, 2]
    min_divE = data[:, 3]
    max_divE = data[:, 4]
    rho_mean = 10 .^ mean_divE
    rho_upper = rho_mean .- 10 .^ (mean_divE .- std_divE)
    rho_lower = 10 .^ (mean_divE .+ std_divE) .- rho_mean

    lines!(ax, 1:length(sims), rho_max, color = 1:length(sims), colormap=:batlow, colorrange=(0,length(sims)+1), linestyle=:dash,
        linewidth=lw_dash)
    band!(ax, 1:length(sims), rho_lower, rho_upper, color=1:length(sims), colormap=:batlow, colorrange=(0,length(sims)+1), alpha=0.8)
    scatter!(ax, 1:length(sims), rho_mean, markersize = 15, strokecolor=:white, strokewidth = 1.0,
            color=1:length(sims), colormap = :batlow, colorrange=(0,length(sims)+1))

    textlabel!(ax, Point2f(length(sims)-0.5, ypos_label), text = "3D",
                    offset = (0, 0), text_color = :black, alpha=0.0, strokewidth=0)

    resize_to_layout!(f)
    save(plot_name, f)
end


function Fig03_1D(data_file, sims, plot_name; y_range=(1.e-12, 1.e-9), ypos_label=3.e-10)

    f = Figure(size=(625, 250))

    lw_dash = 4.0
    yticks = 10.0 .^ (round(Int, log10(y_range[1])):round(Int, log10(y_range[2])))

    ax = Axis(f[1, 1], aspect = 3,
                ylabel=L"\overline{\nabla \cdot \mathbf{E}} \quad [4\pi (\rho_+ - \rho_-)]",
                xticks=1:length(sims),
                xtickformat = x -> sims[Int64.(x)],
                yscale = log10,
                yticks = yticks,
                yminorticks = IntervalsBetween(9),
                xminorticksvisible=false,
                ytickformat = values -> [iszero(value) ? "0" : rich("$(value < 0 ? "-" : "")10", superscript("$(@sprintf("%i", log10(abs(value))))")) for value in values]
                )
    xlims!(ax, (0.5, length(sims) + 0.5))
    ylims!(ax, y_range)

    # columns: mean, std, min, max of log10|divE|
    data = readdlm(data_file)
    mean_divE = data[:, 1]
    std_divE = data[:, 2]
    max_divE = data[:, 4]
    rho_mean = 10 .^ mean_divE
    rho_lower = 10 .^ (mean_divE .- std_divE)
    rho_upper = 10 .^ (mean_divE .+ std_divE)
    rho_max = 10 .^ max_divE

    lines!(ax, 1:length(sims), rho_max, color = 1:length(sims), colormap=:batlow, colorrange=(0,length(sims)+1), linestyle=:dash,
        linewidth=lw_dash)
    band!(ax, 1:length(sims), rho_lower, rho_upper, color=1:length(sims), colormap=:batlow, colorrange=(0,length(sims)+1), alpha=0.8)
    scatter!(ax, 1:length(sims), rho_mean, markersize = 15, strokecolor=:white, strokewidth = 1.0,
            color=1:length(sims), colormap = :batlow, colorrange=(0,length(sims)+1))

    textlabel!(ax, Point2f(length(sims)-0.5, ypos_label), text = "1D",
                    offset = (0, 0), text_color = :black, alpha=0.0, strokewidth=0)

    resize_to_layout!(f)
    save(plot_name, f)
end


set_light_theme!(fontsize=25)

x_labels = [L"\mathcal{Z}", L"\mathcal{S}_1", 
            L"\mathcal{S}_2", L"\mathcal{S}_3", 
            L"\mathcal{S}_4", L"\mathcal{S}_5", 
            L"\mathcal{S}_6", L"\mathcal{S}_7", 
            L"\mathcal{S}_8", L"\mathcal{S}_9",
            L"\mathcal{S}_{10}", L"\mathcal{S}_{11}"
            ]

sim_dims = ["1D", "2D", "3D"]
data_files = [data_path * "charge_conservation/" * sim * "/charge_conservation_error.txt" for sim in sim_dims]

sim_path = data_path * "charge_conservation/"

# for i = 1:length(sim_dims)
#     sim_paths = sim_path * sim_dims[i] * "/" .*
#             ["zigzag", "esirkepov_1st", 
#             "esirkepov_2nd", "esirkepov_3rd", "esirkepov_4th", 
#             "esirkepov_5th", "esirkepov_6th", "esirkepov_7th", "esirkepov_8th",
#             "esirkepov_9th"#, "esirkepov_10th", "esirkepov_11th"
#             ] .* "/quiet_plasma"
#     write_divE_error(sim_paths, data_files[i])
# end

plot_name = plot_path * "Fig03.pdf"
Fig03(data_files, x_labels, plot_name)

#data_file_1D = data_path * "charge_conservation/1D_single/charge_conservation_error.txt"
#plot_name_1D = plot_path * "Fig03.png"
#y_range_1D = (1.e-12, 1.e-9)
#Fig03_1D(data_file_1D, x_labels, plot_name_1D, y_range=y_range_1D)
