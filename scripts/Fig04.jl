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
function Ne2D(ppc)
    cells_per_skindepth = 10.0
    cell_volume = cells_per_skindepth^2
    return cell_volume * ppc/2
end

"""
    3D

The 3D parameter files set ppc0 to half the folder label, so the electron density
per cell is ppc/4 rather than ppc/2.
"""
function Ne3D(ppc)
    cells_per_skindepth = 10.0
    cell_volume = cells_per_skindepth^3
    return cell_volume * ppc/4
end

"""
    Kato 2013, Eq. 34, 2D case
"""
function Neff2D(dEdt)
    return 1.0/(4*abs(dEdt))
end

"""
    Kato 2013, Eq. 34, 3D case
"""
function Neff3D(dEdt, v0; cells_per_skindepth=10.0)
    coulomb_log = log(pi * v0 * cells_per_skindepth)
    return coulomb_log/(4pi * v0 * abs(dEdt))
end

"""
    Speed in units of c of a particle with kinetic energy `E` in units of m_e c^2.
"""
function beta_from_E(E)
    return sqrt(1.0 - 1.0/(1.0 + E)^2)
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
folders = ["zigzag", "esirkepov_2nd", "esirkepov_3rd", "esirkepov_5th", "esirkepov_9th"]
#run_E_data(sim_path, folders)

sim_path = data_path * "stopping_power/3D/"
folders = ["zigzag", "esirkepov_2nd", "esirkepov_3rd"]
#run_E_data(sim_path, folders)


function Fig04(sim_path, sims, titles, dim, plot_name)

    f = Figure(size=(1550, 400))

    folders = ["8ppc", "16ppc", "32ppc", "64ppc"]
    panel_titles = dim == 2 ? ["8ppc", "16ppc", "32ppc", "64ppc"] :
                              ["4ppc", "8ppc", "16ppc", "32ppc"]
    ppc = [8, 16, 32, 64]

    # setup electrons per skin depth volume, and its counterpart inferred from the
    # measured loss rate; the 3D rate carries a velocity and a Coulomb logarithm
    Ne = dim == 2 ? Ne2D : Ne3D
    Neff = dim == 2 ? (dEdt, E0) -> Neff2D(dEdt) :
                      (dEdt, E0) -> Neff3D(dEdt, beta_from_E(E0))

    for i = 1:length(folders)

        println("Processing folder: ", folders[i])

        ax = Axis(f[1, i], xlabel=L"t \quad \left[ \omega_p^{-1} \right]", ylabel=L"\varepsilon \:\: [m_e c^2]",
                    title=L"\text{%$(panel_titles[i])} \rightarrow N_e^{(%$dim)} = %$(round(Ne(ppc[i]); digits=0))",
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

            band!(ax, t, E .- E_std, E .+ E_std, alpha = 0.1, color=j, colormap = :batlow, colorrange=(0,length(sims)+1))

            lines!(ax, t, E,
                color=j, colormap = :batlow, colorrange=(0,length(sims)+1), 
                linewidth=2.5, label=L"N_\mathrm{e,eff} = %$(round( Neff(fit.param[1], E[1]); digits=1))")
        
            lines!(ax, t, E_fit, 
                color=j, colormap = :batlow, colorrange=(0,length(sims)+1), 
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
            
            lins = [lines!(ax, [0.0], [0.0], color=j, colormap = :batlow, colorrange=(0,length(sims)+1), 
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
titles = [L"\text{Zigzag}",
          L"\text{Esirkepov} \:\: \mathcal{S}_2",
          L"\text{Esirkepov} \:\: \mathcal{S}_3",
          L"\text{Esirkepov} \:\: \mathcal{S}_5",
          L"\text{Esirkepov} \:\: \mathcal{S}_9"]
dim = 2
plot_name = plot_path * "Fig04.pdf"
Fig04(sim_path, folders, titles, dim, plot_name)


#sim_path = "/lus/flare/projects/CRsInCosmoSims/PIC/Paper/2600.0000/simulations/stopping_power/3D/"
#folders = ["zigzag", "esirkepov_2nd", "esirkepov_3rd"]
#titles = [L"\text{Zigzag}",
#          L"\text{Esirkepov} \:\: \mathcal{S}_2",
#          L"\text{Esirkepov} \:\: \mathcal{S}_3"]
#dim = 3
#plot_path = "/lus/flare/projects/CRsInCosmoSims/PIC/Paper/2600.0000/Plots/"
#plot_name = plot_path * "Fig04.pdf"
#Fig04(sim_path, folders, titles, dim, plot_name)

