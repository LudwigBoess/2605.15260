include(joinpath(@__DIR__, "config.jl"))

using CairoMakie, MakieHelper
using EntityTools
using LsqFit
using Printf
using LaTeXStrings
using StatsBase
using DelimitedFiles

# Define a model function (e.g., exponential)
model(x, p) = p[1] .* x .+ p[2]


"""
    2D
"""
function Ne(ppc)
    cells_per_skindepth = 10.0
    cell_volume = cells_per_skindepth^2
    return cell_volume * ppc/2
end

function Neff(dEdt)
    return 1.0/(4*abs(dEdt))
end


function write_E_data(sim_path)
    info = EntityData(sim_path)
    E = Vector{Float64}(undef, length(info.steps))
    E_std = Vector{Float64}(undef, length(info.steps))
    for i in 1:length(info.steps)
        p = read_particles(info, 3, ["U1", "U2", "U3"], i_step=info.steps[i])
        _E = sqrt.(p["U1"].^2 .+ p["U2"].^2 .+ p["U3"].^2 .+ 1.0) .- 1.0
        E[i] = mean(_E)
        E_std[i] = std(_E)
    end

    writedlm(sim_path * "/energy_data.txt", [info.times E E_std])
end

function run_E_data(sim_path, sims)


    folders = ["8ppc", "16ppc", "32ppc", 
                "64ppc"]
    for i = 1:length(folders), j = 1:length(sims)

        println("Processing folder: ", sims[j] * "/_$(folders[i])")
        out_path = sim_path * sims[j] * "/_$(folders[i])"

        write_E_data(out_path)
    end

end


sim_path = data_path * "stopping_power/2D/"
folders = ["zigzag", "esirkepov_2nd", "esirkepov_3rd", 
            "esirkepov_5th", 
            "esirkepov_9th"
            ]
#run_E_data(sim_path, folders)


function Fig04(sim_path, sims, plot_name)

    f = Figure(size=(1550, 400))

    titles = [L"\text{Zigzag}", 
            L"\text{Esirkepov} \:\: \mathcal{S}_2", 
            L"\text{Esirkepov} \:\: \mathcal{S}_3", 
            L"\text{Esirkepov} \:\: \mathcal{S}_5",
            L"\text{Esirkepov} \:\: \mathcal{S}_9"]

    folders = ["8ppc", "16ppc", "32ppc", "64ppc"]
    ppc = [8, 16, 32, 64]

    for i = 1:length(folders)

        println("Processing folder: ", folders[i])

        ax = Axis(f[1, i], xlabel=L"t \quad \left[ \omega_p^{-1} \right]", ylabel=L"\varepsilon \:\: [m_e c^2]",
                    title=L"\text{%$(folders[i])} \rightarrow N_e^{(2)} = %$(round(Ne(ppc[i]); digits=0))",
                    aspect = 1)
        xlims!(ax, (0.0, 4.e2))
        #ylims!(ax, (8.2, 9.1))
        ylims!(ax, (8.7, 9.1))
        #ylims!(ax, (7.5, 7.9))

        for j = 1:length(sims)

            out_path = sim_path * sims[j] * "/_$(folders[i])"
            data = readdlm(out_path * "/energy_data.txt")
            t = data[:, 1]
            E = data[:, 2]
            E_std = data[:, 3]

            println("maximum energy: ", maximum(E))
            # Perform the curve fit
            fit = curve_fit(model, t, E, [1.0, 0.0])

            E_fit = Vector{Float64}(undef, length(E))
            for k = 1:length(E)
                E_fit[k] = model(t[k], fit.param)
            end

            band!(ax, t, E .- E_std, E .+ E_std, alpha = 0.1, color=j, colormap = :thermal, colorrange=(0,length(sims)+1))

            lines!(ax, t, E,
                color=j, colormap = :thermal, colorrange=(0,length(sims)+1), 
                linewidth=2.5, label=L"N_\mathrm{e,eff} = %$(round( Neff(fit.param[1]); digits=1))")
        
            lines!(ax, t, E_fit, 
                color=j, colormap = :thermal, colorrange=(0,length(sims)+1), 
                linestyle=:dash, alpha = 0.7,
                linewidth=3.0)
            
            # no energy loss
            hlines!(ax, E[1], color=:black, alpha=0.5, linestyle=:dot)
        end
        axislegend(ax, position = :lb, framevisible = false)

        if i > 1
            hideydecorations!(ax, ticks = false)
        end

        if i == length(folders)
            
            lins = [lines!(ax, [0.0], [0.0], color=j, colormap = :thermal, colorrange=(0,length(sims)+1), 
                    linewidth=2.5, label=titles[j]) for j = 1:length(sims)]
            #f[1, i+1] = Legend(f, ax, framevisible = false)
            Legend(f[1, i+1], lins, titles, framevisible = false)
        end
        data = nothing
    end

    save(plot_name, f)

end

set_light_theme!(fontsize=20)

sim_path = data_path * "stopping_power/2D/"
folders = ["zigzag", "esirkepov_2nd", "esirkepov_3rd", "esirkepov_5th", "esirkepov_9th"]
plot_name = plot_path * "Fig04.pdf"
Fig04(sim_path, folders, plot_name)

