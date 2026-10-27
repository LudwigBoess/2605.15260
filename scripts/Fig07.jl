include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using CairoMakie, MakieHelper
using ColorSchemes
@info "makie done"
using Printf
using EntityTools
using DelimitedFiles
using LinearAlgebra
@info "entitytools done"

# delta/omega_pe
function growth_rate(a, b, c)
    alpha = 1.e-4 # beam to background ratio
    gamma = 100.0 # beam Lorentz factor
    return sqrt(3.0)/2^(4/3) * (alpha/gamma)^(1/3) / sqrt(2)
end

function Fig07(sim_paths, sims, plot_name, t_anchor, E2_anchor, alpha, gamma, L_box)

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
    linestyles = [:dot, :solid]

    for i = 1:length(sim_paths), j = 1:2

        println("Processing folder: ", sims[i])

        d = readdlm(sim_paths[i] * subfolder[j] * "/oblique_modes$(subfolder[j])_stats.csv", ',', header=true)
        t = d[1][:,2]
        E2 = d[1][:,7]

        println(maximum(E2))


        lines!(ax, t, E2, 
            color=i, colormap = :batlow, colorrange=(0,length(sims)-1.5), 
            linewidth=4, linestyle=linestyles[j])
        if j == 2
            lines!(ax, [0.0], [0.0], 
                color=i, colormap = :batlow, colorrange=(0,length(sims)-1.5), 
                linewidth=4, label=sims[i])
        end

        

    end

    # analytic linear growth: field energy grows as exp(2 δ t)
    t_an = LinRange(1200.0, 1700.0, 2)
    la = lines!(ax, t_an, E2_anchor .* exp.(2growth_rate(alpha, gamma, L_box) .* (t_an .- t_anchor)),
                color=:black, linestyle=:dash,
                linewidth=2)

    l0 = lines!(ax, [0.0], [0.0],
                color=:black, linestyle=linestyles[1],
                linewidth=2)
    l8 = lines!(ax, [0.0], [0.0],
                color=:black, linestyle=linestyles[2],
                linewidth=2)

    axislegend(ax; position = :lt, nbanks = 2, framevisible = false)
    axislegend(ax, [l0, l8, la], [L" 0 \: \mathcal{F}", L" 8 \: \mathcal{F}", L"\Gamma_\mathrm{max, obl.}"]; position = :rb, framevisible = false)
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

sim_paths = data_path * "oblique_modes/" .*
            ["zigzag", 
            "esirkepov_2nd", "esirkepov_3rd", 
            "esirkepov_5th", 
            "esirkepov_9th"
            ] .* "/oblique_modes"

plot_name = plot_path * "Fig07.pdf"
t_anchor  = 1200.0   # [ω_p^-1] anchor of the analytic growth line
E2_anchor = 3.e-9    # [n_b m_e c^2] E_y^2/8π at t_anchor
alpha     = 1.e-4                # beam to background density ratio (pair totals)
gamma     = sqrt(1 + 100.0^2)    # beam Lorentz factor, drift 4-velocity u = 100
L_box     = 25.6                 # [c/ω_p] box length along the beam
Fig07(sim_paths, sims, plot_name, t_anchor, E2_anchor, alpha, gamma, L_box)

