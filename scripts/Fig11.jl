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



function Fig11(sim_paths, sims, plot_name)

    f = Figure(size=(585, 1235))

    ax1 = Axis(f[2, 1:6], xlabel=L"\text{Time} \quad \omega_e t", ylabel=L"\varepsilon \:\: [m_e c^2]",
                aspect = 1,
                xticks = [0.0, 250.0, 500.0, 750.0, 1000.0]
                )
    xlims!(ax1, (0.0, 1100.0))
    ylims!(ax1, (-2.0, 4.0))
    hidexdecorations!(ax1, ticks = false, minorticks = false)

    ax_zoom_l = Axis(f[1, 1:3], yticklabelsize=15,
                    yticks = [2.95, 3.0, 3.05, 3.1, 3.15])
    xlims!(ax_zoom_l, (0, 100))
    ylims!(ax_zoom_l, (3.0, 3.10))
    zoom_lines!(ax1, ax_zoom_l,
                rectattrs=(linestyle=:dash, strokewidth=1.0, transparency=0.5),
                lineattrs=(linestyle=:dash, linewidth=1.5, transparency=0.5))
    hidexdecorations!(ax_zoom_l, ticks = false, minorticks = false)
    hideydecorations!(ax_zoom_l, ticks = false, minorticks = false)
                

    ax_zoom_r = Axis(f[1, 4:6], yticklabelsize=15, yticklabelpad=15,
                yticks = [2.95, 3.0, 3.05, 3.1, 3.15])
    xlims!(ax_zoom_r, (1000, 1100))
    ylims!(ax_zoom_r, (3.0, 3.1))
    zoom_lines!(ax1, ax_zoom_r,
                rectattrs=(linestyle=:dash, strokewidth=1.5, transparency=0.5),
                lineattrs=(linestyle=:dash, linewidth=1.5, transparency=0.5))
    hidexdecorations!(ax_zoom_r, ticks = false, minorticks = false)


    ax2 = Axis(f[3, 1:6], xlabel=L"\text{Time} \quad t \quad \left[ \omega_p^{-1} \right]", 
            ylabel=L"\log_{10} \:\: | \varepsilon -  \varepsilon_0 | / \varepsilon_0",
                aspect = 1.2,
                xticks = [0.0, 250.0, 500.0, 750.0, 1000.0]
                )
    xlims!(ax2, (0.0, 1100.0))
    ylims!(ax2, (-5, -1))

    rowsize!(f.layout, 1, 150.0)
    rowgap!(f.layout, 1, 10.0)

    E0 = 0.0

    for i = 1:5

        println("Processing folder: ", sims[i])

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

        if i == 1
            E0 = Etot[1]
        end

        lines!(ax1, t, Etot, 
            color=i, colormap = :batlow, colorrange=(0,6), 
            linewidth=2.5, label=sims[i])

        lines!(ax_zoom_l, t, Etot, 
            color=i, colormap = :batlow, colorrange=(0,6), 
            linewidth=2.5, label=sims[i])

        lines!(ax_zoom_r, t, Etot, 
            color=i, colormap = :batlow, colorrange=(0,6), 
            linewidth=2.5, label=sims[i])

        if i == 1
            lines!(ax_zoom_l, [0.0, t[end]], [Etot[1], Etot[1]],  
                color=:black, linestyle=:dot, alpha=0.7,
                linewidth=2.5)

            lines!(ax_zoom_r, [0.0, t[end]], [Etot[1], Etot[1]],  
                color=:black, linestyle=:dot, alpha=0.7,
                linewidth=2.5)
        end

        lines!(ax1, t, B, 
            color=i, colormap = :batlow, colorrange=(0,6), 
            linewidth=2.5, linestyle=:dash)

        lines!(ax1, t, E, 
            color=i, colormap = :batlow, colorrange=(0,6), 
            linewidth=2.5, linestyle=:dot)

        lines!(ax1, t, T, 
            color=i, colormap = :batlow, colorrange=(0,6), 
            linewidth=2.5, linestyle=:dashdot)


        E_err = log10.(abs.(Etot .- Etot[1]) ./ Etot[1])
        lines!(ax2, t, E_err, 
            color=i, colormap = :batlow, colorrange=(0,6), 
            linewidth=2.5)

        println("Max energy error: ", maximum(E_err))
            
        data = nothing
    end

    linestyles = [:solid, :dash, :dot, :dashdot, :dashdotdot]

    for i = 6:length(sim_paths)
        println("Processing folder: ", sims[i])

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

        p = sortperm(t)

        Etot = 0.5 .* (B .+ E) .+ T
        E_err = log10.(abs.(Etot .- Etot[1]) ./ Etot[1])
        lines!(ax2, t[p], E_err[p], 
            color=1, colormap = :batlow, colorrange=(0,6), 
            linewidth=2.5, linestyle=linestyles[i-4])
    end

    legend_labels = [L"\Delta x_0",
        L"\Delta x_0 +16 \mathcal{F}",
        L"\frac{1}{2} \Delta x_0",
        L"\frac{1}{4} \Delta x_0",
        L"\frac{1}{8} \Delta x_0"]

    for i = 5:length(sim_paths)
        lines!(ax2, [0.0], [0.0], 
            color=:black, 
            linewidth=2, linestyle=linestyles[i-4],
            label=legend_labels[i-4])
    end

    hlines!(ax_zoom_l, E0, 
            color=:black, linewidth=2.0, linestyle=:dash, 
            alpha=0.5)

    hlines!(ax_zoom_r, E0, 
            color=:black, linewidth=2.0, linestyle=:dash, 
            alpha=0.5)

    lines!(ax1, [0.0], [0.0], 
            color=:black,
            linewidth=2.5, label=L"\varepsilon_\text{tot}")
    lines!(ax1, [0.0], [0.0], 
        color=:black,
        linewidth=2.5, linestyle=:dash, label=L"\varepsilon_\text{B}")
    lines!(ax1, [0.0], [0.0], 
        color=:black,
        linewidth=2.5, linestyle=:dot, label=L"\varepsilon_\text{E}")
    lines!(ax1, [0.0], [0.0], 
        color=:black,
        linewidth=2.5, linestyle=:dashdot, label=L"\varepsilon_\text{kin}")

    axislegend(ax1, position = :cb, nbanks=5, framevisible = false, 
            labelsize=30)
    axislegend(ax2, position = :lt, nbanks=5, framevisible = false, 
            labelsize=20)
    resize_to_layout!(f)
    save(plot_name, f)

end

set_light_theme!(fontsize=30)

runs = [L"\mathcal{Z}", 
        L"\mathcal{S}_2", 
        L"\mathcal{S}_3", 
        L"\mathcal{S}_5",
        L"\mathcal{S}_9",
        L"\Delta x_0",
        L"\Delta x_0 +16 \mathcal{F}",
        L"\frac{1}{2} \Delta x_0",
        L"\frac{1}{4} \Delta x_0",
        L"\frac{1}{8} \Delta x_0"
        ]

sim_paths = data_path * "decaying_turb/2D/" .*
            ["zigzag", "esirkepov_2nd", "esirkepov_3rd", "esirkepov_5th", "esirkepov_9th", 
            "zigzag_16filters", "zigzag_2x", "zigzag_4x", "zigzag_8x"
            ] .* "/turbulence"


t = 1092.0
plot_name = plot_path * "Fig11.pdf"
Fig11(sim_paths, runs, plot_name)
