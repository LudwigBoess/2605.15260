include(joinpath(@__DIR__, "config.jl"))

using CairoMakie, MakieHelper
using Printf
using DelimitedFiles


function Fig09(sim_paths, sims, plot_name)

    f = Figure(size=(600, 500))

    ax = Axis(f[1, 1], xlabel=L"t \quad \left[10^4 \: \omega_p^{-1}\right]", ylabel=L"\log_{10} \langle B_{yz} \rangle \:\: \left[ \sqrt{4\pi n_e m_e c^2} \right]",
                aspect = 1.2,
                )
    xlims!(ax, (0.0, 15.0))
    ylims!(ax, (-7, 2))

    for i = 1:length(sim_paths)

        println("Processing folder: ", sims[i])

        out_path = sim_paths[i]
        if !isdir(out_path)
            continue
        end
        fi = out_path * "/_bell_stats.csv"
        d = readdlm(fi, ',', header=true)
        t = d[1][:,2] .* 1.e-4
        By = d[1][:,4]
        Bz = d[1][:,5]

        Bxy = sqrt.(By.^2 .+ Bz.^2)

        maximum(Bxy) |> println

        if i < length(sim_paths)
            lines!(ax, t, log10.(Bxy), 
            color=i, colormap = :thermal, colorrange=(0,length(sims)), 
            linewidth=2)
            lines!(ax, [0.0], [0.0], 
            color=i, colormap = :thermal, colorrange=(0,length(sims)), 
            linewidth=4, label=sims[i])
        else
            lines!(ax, t, log10.(Bxy), 
            color=:black, linestyle=:dot, 
            linewidth=2)
            lines!(ax, [0.0], [0.0], 
            color=:black, linestyle=:dot, 
            linewidth=4, label=sims[i])
        end
        

    end
    axislegend(ax; position = :rb, nbanks = 2, framevisible = false)
    resize_to_layout!(f)
    save(plot_name, f)
end

set_light_theme!(fontsize=25)

sims = ["Zigzag", 
            L"\text{Esirkepov} \:\: \mathcal{S}_1", 
            L"\text{Esirkepov} \:\: \mathcal{S}_2", 
            L"\text{Esirkepov} \:\: \mathcal{S}_3", 
            L"\text{Esirkepov} \:\: \mathcal{S}_4", 
            L"\text{Esirkepov} \:\: \mathcal{S}_5",
            L"\text{Esirkepov} \:\: \mathcal{S}_6", 
            L"\text{Esirkepov} \:\: \mathcal{S}_7", 
            L"\text{Esirkepov} \:\: \mathcal{S}_8", 
            L"\text{Esirkepov} \:\: \mathcal{S}_9",
            L"\text{Zigzag} + 8 \mathcal{F}",
            ]


sim_paths = data_path * "bell/1D/" .*
            ["zigzag", "esirkepov_1st", 
            "esirkepov_2nd", "esirkepov_3rd", "esirkepov_4th", 
            "esirkepov_5th", "esirkepov_6th", "esirkepov_7th", "esirkepov_8th",
            "esirkepov_9th",
            "zigzag_8filters",
            ] .* "/_bell"

plot_name = plot_path * "Fig09.pdf"
Fig09(sim_paths, sims, plot_name)