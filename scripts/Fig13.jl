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


function write_powerspec(sim_paths, t, cells_per_de)

    for i = 1:length(sim_paths)
        @info "reading data 2D"
        info = EntityData(sim_paths[i])
        # power spec
        d = read_field(info, ["B1", "B2", "B3"], verbose=true; t)

        k, Pk = power_spectrum(d["B1"], d["B2"], d["B3"], cells_per_skindepth = cells_per_de[i])

        writedlm(sim_paths[i] * "/powerspec.txt", [k Pk])
    end
end

function add_analytic!(ax)
    # Define k0 and Pk0 for the k^(-5/3) reference line
    # Using the x-axis limits for the start and end points of the line
    k0 = [4e-3, 2.0e-2]
    
    # C is an arbitrary constant to vertically shift the line next to your data.
    # You may need to tweak this value (e.g., 1e6, 1e8) depending on your spectrum's amplitude.
    C = 1.0e-4
    Pk0 = C .* k0 .^ (-5.0 / 3.0) .* k0[2]

    k1 = [5.e-2, 0.2]
    C2 = 1.0e-8
    Pk1 = C2 .* k1 .^ (-4) .* k1[2]


    lines!(ax, k0, Pk0, color=:black, linestyle=:dash, linewidth=2.0)
    lines!(ax, k1, Pk1, color=:black, linestyle=:dot, linewidth=2.0)
    text!(ax, 3.0e-3, C * 10 * 1.0^(-5.0/3.0), text = L"k^{-5/3}", color=:black, align = (:left, :bottom), offset=(5, -5))
    text!(ax, 3.0e-2, C2 * 600 * 1.0^(-4), text = L"k^{-4}", color=:black, align = (:left, :bottom), offset=(5, -5))
end

function Fig13(sim_paths, run_labels, t, factors, plot_name)

    f = Figure(size=(550, 1200))
    box_with = 512.0

    # top panel
    ax1 = Axis(f[1, 1], xlabel = L"k_\perp d_e", ylabel=L"B^2(k_\perp)\:\: [a.u.]", aspect = 1.3,
        xscale=log10, yscale=log10, 
        xminorticks = IntervalsBetween(9),
        yminorticks = IntervalsBetween(9))
    xlims!(ax1, (1.e-3, 100.0) )
    ylims!(ax1, (3.e-7,  1.0) )
    hidexdecorations!(ax1, ticks = false, minorticks=false)

    for i = 1:length(sim_paths)
        d = readdlm(sim_paths[i] * "/turbulence/powerspec.txt")
        k = d[:, 1]
        Pk = d[:, 2] ./ d[1,2]

        lines!(ax1, k, Pk, color=i, colormap = :thermal, colorrange=(0,length(sim_paths)+1), 
            linewidth=2.5, label=run_labels[i])
    end

    zigzag_runs = ["_16filters", "_2x", "_4x", "_8x"]
    subrun_labels = [L"\mathcal{Z}; 16 \mathcal{F}", L"\mathcal{Z}; \frac{1}{2}\Delta x_0", L"\mathcal{Z}; \frac{1}{4}\Delta x_0", L"\mathcal{Z}; \frac{1}{8}\Delta x_0"]
    for i = 1:length(zigzag_runs)
        d = readdlm(sim_paths[1] * zigzag_runs[i] * "/turbulence/powerspec.txt")
        k = d[:, 1]
        Pk = d[:, 2] ./ d[1,2]


        println("max Pk: ", maximum(Pk))
        println("min k: ", minimum(k), " max k: ", maximum(k))

        lines!(ax1, k, Pk, color=i, colormap = :viridis, colorrange=(0,4), 
            linewidth=2.5, label=subrun_labels[i])
    end
    add_analytic!(ax1)
    axislegend(ax1, framevisible=false, position = :rt, nbanks = 2, labelsize=20)

    # middle panel
    ax2 = Axis(f[2,1], xlabel = L"k_\perp d_e", ylabel=L"B^2(k_\perp)\:\: [a.u.]", aspect = 1.3,
        xscale=log10, yscale=log10, 
        xminorticks = IntervalsBetween(9),
        yminorticks = IntervalsBetween(9))
    xlims!(ax2, (1.e-3, 100.0) )
    ylims!(ax2, (3.e-7,  1.0) )
    hidexdecorations!(ax2, ticks = false, minorticks=false)

    zigzag_runs = ["", "_1filters", "_2filters", "_4filters", "_8filters", "_16filters"]
    subrun_labels = [L"0 \mathcal{F}", L"1 \mathcal{F}", L"2 \mathcal{F}", L"4 \mathcal{F}", L"8 \mathcal{F}", L"16 \mathcal{F}"]
    for i = 1:length(zigzag_runs)
        d = readdlm(sim_paths[1] * zigzag_runs[i] * "/turbulence/powerspec.txt")
        k = d[:, 1]
        Pk = d[:, 2] ./ d[1,2]

        lines!(ax2, k, Pk, color=i, colormap = :thermal, colorrange=(0,7), 
            linewidth=2.5, label=subrun_labels[i])
    end
    add_analytic!(ax2)
    axislegend(ax2, framevisible=false, position = :rt, nbanks = 2)
    text!(ax2, 0.7, 0.3, text = "Zigzag", color=:black, align = (:right, :center), offset=(5, -5))

    # bottom panel
    ax3 = Axis(f[3, 1], xlabel = L"k_\perp d_e", ylabel=L"B^2(k_\perp)\:\: [a.u.]", aspect = 1.3,
        xscale=log10, yscale=log10, 
        xminorticks = IntervalsBetween(9),
        yminorticks = IntervalsBetween(9))
    xlims!(ax3, (1.e-3, 100.0) )
    ylims!(ax3, (3.e-7,  1.0) )

    zigzag_runs = ["", "_1filters", "_2filters", "_4filters", "_8filters", "_16filters"]
    subrun_labels = [L"0 \mathcal{F}", L"1 \mathcal{F}", L"2 \mathcal{F}", L"4 \mathcal{F}", L"8 \mathcal{F}", L"16 \mathcal{F}"]
    for i = 1:length(zigzag_runs)
        d = readdlm(sim_paths[2] * zigzag_runs[i] * "/turbulence/powerspec.txt")
        k = d[:, 1]
        Pk = d[:, 2] ./ d[1,2]

        lines!(ax3, k, Pk, color=i, colormap = :thermal, colorrange=(0,7), 
            linewidth=2.5, label=subrun_labels[i])
    end
    add_analytic!(ax3)
    axislegend(ax3, framevisible=false, position = :rt, nbanks = 2)
    text!(ax3, 1.0, 0.3, text = L"\text{Esirkepov} \: \mathcal{S}_2", color=:black, align = (:right, :center), offset=(5, -5))

    resize_to_layout!(f)
    @info "saving"
    save(plot_name, f)

end

set_light_theme!(fontsize=24)


runs = [L"\mathcal{Z}", 
        L"\mathcal{S}_2", 
        L"\mathcal{S}_3", 
        L"\mathcal{S}_5",
        L"\mathcal{S}_9"]

sim_paths = data_path * "decaying_turb/2D/" .*
            ["zigzag", "esirkepov_2nd", "esirkepov_3rd", "esirkepov_5th", "esirkepov_9th"]

factors = 4 .* [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 4.0, 8.0]
t = 1092.0
write_powerspec(sim_paths, t, factors)

plot_name = plot_path * "Fig13.pdf"

Fig13(sim_paths, runs, t, factors, plot_name)