include(joinpath(@__DIR__, "config.jl"))

using CairoMakie, MakieHelper
using Printf
using DelimitedFiles


# Bell (non-resonant) cosmic-ray current-driven instability (Bell 2004,
# MNRAS 353, 550). In entity code units (c = 1, larmor0 = skindepth0 = 1
# ⇒ B0 = 1, σ_0 = 1, and the 4π/c of Ampère's law absorbed into J) the
# dispersion relation is
#
#     γ(k)² = k v_A² (ζ − k),     ζ = J_cr / B0
#
# so the fastest-growing mode sits at k_max = ζ/2 with γ_max = v_A k_max, and
# modes above k_c = ζ are stable. Rates come out in units of the fiducial plasma
# frequency ω_p0 — the unit of the "time" column of the *_stats.csv files.
function bell_growth_rate(; density_cr = 0.01,
                            u_cr       = 99.9962499609,
                            density_bg = 1.0,
                            m_i        = 50.0,
                            m_e        = 1.0,
                            B0         = 1.0)
    n_cr = 0.5 * density_cr
    n_bg = 0.5 * density_bg
    v_cr = u_cr / sqrt(1 + u_cr^2)
    J_cr = n_cr * v_cr
    v_A  = B0 / sqrt(n_bg * (m_i + m_e))
    return 0.5 * v_A * J_cr / B0
end


function Fig09(sim_paths, sims, plot_name)

    f = Figure(size=(600, 500))

    ax = Axis(f[1, 1], xlabel=L"t \quad \left[10^4 \: \omega_p^{-1}\right]", ylabel=L"\log_{10} \langle B_{yz} \rangle \:\: \left[ \sqrt{4\pi n_e m_e c^2} \right]",
                aspect = 1.2,
                )
    xlims!(ax, (0.0, 15.0))
    ylims!(ax, (-4, 1))
    #ylims!(ax, (-5, 1))

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

        Bxy = sqrt.(By .+ Bz)

        maximum(Bxy) |> println

        if i < length(sim_paths)
            lines!(ax, t, log10.(Bxy), 
            color=i, colormap = :batlow, colorrange=(0,length(sims)), 
            linewidth=2)
            lines!(ax, [0.0], [0.0], 
            color=i, colormap = :batlow, colorrange=(0,length(sims)), 
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

    # ideal-MHD maximum rate, anchored to the start of the fitted line so the two
    # slopes can be read off against each other
    γ_max    = bell_growth_rate()
    t_ref    = [1.0, 2.0] 
    logB_ref = -2.0
    logB_max = logB_ref[1] .+ (γ_max / log(10)) .* (t_ref .- t_ref[1]) .* 1.e4
    lines!(ax, t_ref, logB_max, color=:black, linestyle=:dash, linewidth=2.5,
            label=L"\Gamma_\mathrm{max,Bell}")

    axislegend(ax; position = :rb, nbanks = 2, framevisible = false)
    resize_to_layout!(f)
    save(plot_name, f)
end



set_light_theme!(fontsize=25)

sims = ["Zigzag", 
            L"\text{Esirkepov} \:\: \mathcal{S}_1", 
            L"\text{Esirkepov} \:\: \mathcal{S}_2", 
            L"\text{Esirkepov} \:\: \mathcal{S}_3", 
            L"\text{Esirkepov} \:\: \mathcal{S}_5",
            L"\text{Esirkepov} \:\: \mathcal{S}_9",
            L"\text{Zigzag} + 8 \mathcal{F}",
            ]


sim_paths = data_path * "bell/" .*
            ["zigzag", 
            "esirkepov_1st", 
            "esirkepov_2nd", 
            "esirkepov_3rd",
            "esirkepov_5th",
            "esirkepov_9th",
            "zigzag_8filters",
            ] .* "/_bell"

plot_name = plot_path * "Fig09.pdf"
Fig09(sim_paths, sims, plot_name)