include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using CairoMakie, MakieHelper
using ColorSchemes
@info "makie done"
using Printf
using EntityTools
using DelimitedFiles
@info "entitytools done"


function Fig07(sim_paths, sims, plot_name)

    f = Figure(size=(650, 600), figure_padding=(5, 30, 5, 15))

    ax = Axis(f[1, 1], xlabel=L"\text{Time} \:\: t \:\: \left[ \omega_p^{-1} \right]", ylabel=L"\langle E_y^2 / 8\pi \rangle \:\: \left[ n_b m_e c^2 \right]",
                aspect = 1,
                yscale = log10,
                yminorticks = IntervalsBetween(9),
                yticks = (10.0.^LinRange(-10, 0, 11), [L"10^{-10}", L"10^{-9}", L"10^{-8}", L"10^{-7}", L"10^{-6}", L"10^{-5}", L"10^{-4}", L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^{0}"])
                )
    xlims!(ax, (0.0, 3000.0))
    ylims!(ax, (1.e-10, 1.e-2))
    subfolder = ["", "_8filters"]
    linestyles = [:dash, :solid]

    for i = 1:length(sim_paths), j = 1:2

        println("Processing folder: ", sims[i])

        d = readdlm(sim_paths[i] * subfolder[j] * "/oblique_modes$(subfolder[j])_stats.csv", ',', header=true)
        t = d[1][:,2]
        E2 = d[1][:,7]

        println(maximum(E2))


        lines!(ax, t, E2, 
            color=i, colormap = :thermal, colorrange=(0,length(sims)-1.5), 
            linewidth=4, linestyle=linestyles[j])
        if j == 2
            lines!(ax, [0.0], [0.0], 
                color=i, colormap = :thermal, colorrange=(0,length(sims)-1.5), 
                linewidth=4, label=sims[i])
        end

        

    end

    l0 = lines!(ax, [0.0], [0.0],
                color=:black, linestyle=linestyles[1],
                linewidth=2)
    l8 = lines!(ax, [0.0], [0.0],
                color=:black, linestyle=linestyles[2],
                linewidth=2)

    axislegend(ax; position = :lt, nbanks = 2, framevisible = false)
    axislegend(ax, [l0, l8], [L" 0 \: \mathcal{F}", L" 8 \: \mathcal{F}"]; position = :rb, framevisible = false)
    resize_to_layout!(f)
    save(plot_name, f)
end

set_light_theme!(fontsize=27)

sims = ["Zigzag", 
            L"\text{Esirkepov} \:\: \mathcal{S}_2", 
            L"\text{Esirkepov} \:\: \mathcal{S}_3", 
            L"\text{Esirkepov} \:\: \mathcal{S}_5", 
            L"\text{Esirkepov} \:\: \mathcal{S}_9",
            L"\text{Zigzag} + 8 \mathcal{F}",
            L"\text{Esirkepov} \:\: \mathcal{S}_{2} + 8 \mathcal{F}"
            ]


sim_paths = data_path * "oblique_modes_backup/" .*
            ["zigzag", 
            "esirkepov_2nd", "esirkepov_3rd", 
            "esirkepov_5th", 
            "esirkepov_9th"
            ] .* "/oblique_modes"

plot_name = plot_path * "Fig07.pdf"
Fig07(sim_paths, sims, plot_name)

