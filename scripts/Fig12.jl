include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using CairoMakie,MakieHelper
using MakieExtra
using ColorSchemes, Colors
@info "makie done"
using Printf
using EntityTools
using DelimitedFiles
@info "entitytools done"

function Fig12(sim_paths, runs, run_colors, plot_name)

    myscale = Symlog10(1.e-4, linscale=0.5)

    f = Figure(size=(800, 600))

    ax1 = Axis(f[1, 1],  xlabel=L"\text{Filter passes} \:\: \mathcal{F}", ylabel=L"(\varepsilon -  \varepsilon_0) / \varepsilon_0",
                aspect = 1.3,
                yminorticks=IntervalsBetween(9),
                xticks=1:length(runs),
                xtickformat = x -> runs[Int64.(x)],
                xminorticksvisible=false,
                yscale=myscale,
                yticks = [-1.0e-1, -1.0e-2, -1.0e-3, -1.0e-4, 0.0, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1],
                ytickformat = values -> [iszero(value) ? "0" : rich("$(value < 0 ? "-" : "")10", superscript("$(@sprintf("%i", log10(abs(value))))")) for value in values]
                )
    xlims!(ax1, (0.5, length(runs)+0.5))
    ylims!(ax1, (-1.0e-1, 1.0e-1))

        # zigzag runs
    band!(ax1, [0.5, 6.5], -1, 1, 
        color=1, colormap=:thermal, colorrange=(0,13), 
        alpha=0.3)

    # 4
    textlabel!(ax1, Point2f(3.5, 3e-2), text = L"\mathcal{Z}",
                    offset = (0, 0), text_color = :black, 
                    alpha=0.75, strokewidth=0)

    band!(ax1, [6.5, 12.5], -1, 1, 
        color=3, colormap=:thermal, colorrange=(0,13), 
        alpha=0.3)

    textlabel!(ax1, Point2f(9, 3e-2), text = L"\mathcal{S}_2",
                    offset = (0, 0), text_color = :black, 
                    alpha=0.75, strokewidth=0)

    band!(ax1, [12.5, 13.5], -1, 1, 
        color=4, colormap=:thermal, colorrange=(0,13), 
        alpha=0.3)
    textlabel!(ax1, Point2f(13, 3e-2), text = L"\mathcal{S}_3",
                    offset = (0, 0), text_color = :black, 
                    alpha=0.75, strokewidth=0)

    band!(ax1, [13.5, 14.5], -1, 1, 
        color=6, colormap=:thermal, colorrange=(0,13), 
        alpha=0.3)
    textlabel!(ax1, Point2f(14, 3e-2), text = L"\mathcal{S}_5",
                    offset = (0, 0), text_color = :black, 
                    alpha=0.75, strokewidth=0)

    band!(ax1, [14.5, 15.5], -1, 1, 
        color=10, colormap=:thermal, colorrange=(0,13), 
        alpha=0.3)
    textlabel!(ax1, Point2f(15, 3e-2), text = L"\mathcal{S}_9",
                    offset = (0, 0), text_color = :black, 
                    alpha=0.75, strokewidth=0)


    for i = 1:length(sim_paths)

        println("Processing folder: ", sim_paths[i])

        out_path = sim_paths[i]
        if !isdir(out_path)
            continue
        end
        fi = out_path * "/turbulence_stats.csv"
        d = readdlm(fi, ',', header=true)
        t = d[1][:,2]
        B = d[1][:,3] + d[1][:,4] + d[1][:,5]
        E = d[1][:,6] + d[1][:,7] + d[1][:,8]
        T = d[1][:,end-1]

        Etot = 0.5 .* (B .+ E) .+ T

        E_err = (Etot[end] - Etot[1]) ./ Etot[1]

        println("Energy error for run: ", E_err)

        scatter!(ax1, i, E_err, 
            markersize = 30, strokecolor=:white, strokewidth = 1.0,
            color=run_colors[i], colormap = :thermal, colorrange=(0,13))

        scatter!(ax1, i, abs.(E_err),
            markersize = 30, 
            strokecolor=:white, strokewidth = 2.0,
            color=:transparent)
    end

    resize_to_layout!(f)
    save(plot_name, f)
end


runs = [L"\mathcal{Z}", 
        L"1 \mathcal{F}",
        L"2 \mathcal{F}",
        L"4 \mathcal{F}",
        L"8 \mathcal{F}",
        L"16 \mathcal{F}",
        L"\mathcal{S}_2", 
        L"1 \mathcal{F}",
        L"2 \mathcal{F}",
        L"4 \mathcal{F}",
        L"8 \mathcal{F}",
        L"16 \mathcal{F}",
        L"\mathcal{S}_3", 
        L"\mathcal{S}_5",
        L"\mathcal{S}_9"
        ]

runs = [L"0", 
        L"1",
        L"2",
        L"4",
        L"8",
        L"16",
        L"0", 
        L"1",
        L"2",
        L"4",
        L"8",
        L"16",
        L"0", 
        L"0",
        L"0"
        ]

run_colors = [1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 4, 6, 10]

sim_paths = data_path * "decaying_turb/2D/" .*
            ["zigzag", "zigzag_1filters", "zigzag_2filters", "zigzag_4filters", "zigzag_8filters", "zigzag_16filters",
            "esirkepov_2nd", "esirkepov_2nd_1filters", "esirkepov_2nd_2filters", "esirkepov_2nd_4filters", "esirkepov_2nd_8filters", "esirkepov_2nd_16filters",
            "esirkepov_3rd", "esirkepov_5th", "esirkepov_9th"
            ] .* "/turbulence"

plot_name = plot_path * "Fig12.pdf"
Fig12(sim_paths, runs, run_colors, plot_name)