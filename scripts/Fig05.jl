include(joinpath(@__DIR__, "config.jl"))

using CairoMakie, MakieHelper
using EntityTools
using StatsBase
using Printf
using ProgressMeter
using DelimitedFiles

function Fig05(sim_paths, time, plot_name)

    t_lim = (0.0, 10000.0)

    f = Figure(size=(1450, 380))
    
    cells = ["0-5", "1", "2", "4", "8", "16"] .* "cell"

    titles = [L"\text{Zigzag}", 
            L"\text{Esirkepov} \:\: \mathcal{S}_2", 
            L"\text{Esirkepov} \:\: \mathcal{S}_3", 
            L"\text{Esirkepov} \:\: \mathcal{S}_5",
            L"\text{Esirkepov} \:\: \mathcal{S}_{9}"]

    labels = [L"d_e = \frac{1}{2} \: \Delta x ",
            L"d_e = \Delta x ", 
            L"d_e = 2 \: \Delta x", 
            L"d_e = 4 \: \Delta x", 
            L"d_e = 8 \: \Delta x",
            L"d_e = 16 \: \Delta x",
            L"0 \: \mathcal{F}", L"8 \: \mathcal{F}"]

    filters = ["", "_8filter"]
    linestyles = [:dash, :solid]
    lw = 3

    for i = 1:length(sim_paths)
        println("reading data from $(sim_paths[i])")

        ax1 = Axis(f[1,i], xlabel = L"t \:\: [10^3 \: \omega_p^{-1}]", ylabel = L"\log_{10}((T - T_0)/T_0)",
                    title=titles[i],        
                    yticks = -6:5, xticks = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0],
                aspect=0.8)

        xlims!(ax1, (0, 10.0))
        ylims!(ax1, (-6, 5))

        for j = 1:2
            for k = 1:length(cells)
                
                #fi = sim_paths[i] * filters[j] * "/heating_" * cells[k]
                #info = EntityData(fi)
                #data = read_field(info, ["T22"], t=[0.0, Inf], reduction_function=sum, verbose=false)
                #println(maximum(data["T22"]))

                #writedlm(sim_paths[i] * filters[j] * "/heating_" * cells[k] * "_" * filters[j] * "_T22.txt", hcat(data["t"], data["T22"]))
                d = readdlm(sim_paths[i] * filters[j] * "/heating_" * cells[k] * "_" * filters[j] * "_T22.txt")
                data = Dict("t" => d[:,1], "T22" => d[:,2])

                rel_err = log10.(abs.(data["T22"] .- data["T22"][1]) ./ data["T22"][1])
                println(maximum(rel_err))
                lines!(ax1, 1.e-3 .* data["t"], rel_err, color = k, colormap=:thermal, colorrange=(0,length(cells)+1), linewidth=lw,
                        linestyle=linestyles[j], label = labels[k])

                data = nothing 
                GC.gc()
            end
        end

        if i > 1
            hideydecorations!(ax1, ticks = false, minorticks=false)        
            if i == length(sim_paths)
                lins = vcat([lines!(ax1, [0.0], [0.0], color = k, colormap=:thermal, colorrange=(0,length(cells)+1), linewidth=lw) 
                                for k = 1:length(cells)],
                            [lines!(ax1, [0.0], [0.0], color = :black, linewidth=2.0, linestyle=linestyles[j]) 
                                for j = 1:2])
                
                Legend(f[1, length(cells)], lins, labels, framevisible = false, orientation = :horizontal, 
                        nbanks=length(cells) + 2, tellheight = true)
            end
        end
    end
    
    save(plot_name, f)

end

set_light_theme!(fontsize=22)

sim_paths = data_path * "numerical_heating/2D/" .*
            ["zigzag", 
            "esirkepov_2nd", 
            "esirkepov_3rd", 
            "esirkepov_5th", 
            "esirkepov_9th"
                ]
plot_name = plot_path * "Fig05.pdf"
Fig05(sim_paths, 10000, plot_name)
