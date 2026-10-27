include(joinpath(@__DIR__, "config.jl"))

using CairoMakie, MakieHelper
using Printf
using StatsBase
using EntityTools
using DelimitedFiles


function read_E(info, t, Nspecies, xmin, xmax)

    data = read_particles(info, Nspecies, ["X1", "U1", "U2", "U3"], t=t, #slice=[[xmin, xmax], [-Inf, Inf], [-Inf, Inf]], 
                                                                    verbose=true)

    sel = findall(xmin .< data["X1"] .< xmax)
    γ = @. √(1 + data["U1"][sel]^2 + data["U2"][sel]^2 + data["U3"][sel]^2) - 1.0

    return γ
end

"""
    get_histograms(bins, Q)

Computes the histograms for `M_s` of all shocked particles
"""
function get_histograms(bins, Q)

    # fit histograms
    hist = fit(Histogram, Q, bins)

    return hist.weights ./ sum(hist.weights)
end


"""
    construct_bins_and_centers(bin_min, bin_max, Nbins)

Returns the bins and their center points.
"""
function construct_bins_and_centers(bin_min, bin_max, Nbins)

    # construct bin boundaries
    bins = LinRange(bin_min, bin_max, Nbins + 1)

    # construct bin centers
    bin_centers = Vector{Float64}(undef, Nbins)

    for i = 1:Nbins
        bin_centers[i] = 0.5 * (bins[i] + bins[i+1])
    end

    bins, bin_centers
end

function write_spectra(sim_path, time)

    shapes = ["zigzag", "esirkepov_3rd"]
    fieldsolvers = ["_Yee", "_min3"]

    downstr_xmin = 1000.0
    downstr_xmax = 2000.0

    upstr_xmin = 4000.0
    upstr_xmax = 5000.0
    Nbins = 900

    bins, bin_centers = construct_bins_and_centers(log10(1.e0), log10(1.e3), Nbins)

    bin_centers = 10.0 .^ bin_centers
    bins = 10.0 .^ bins

    dγ = [bins[i+1] - bins[i] for i = 1:length(bins)-1]

    sim = joinpath(sim_path, shapes[2] * "/" * fieldsolvers[2])
    info = EntityData(sim)

    γm1 = read_E(info, 0.0, 1, upstr_xmin, upstr_xmax)
    γm1_hist = get_histograms(bins, γm1)
    spectrum = (bin_centers).^2 .* γm1_hist ./ dγ

    writedlm(sim_path * "/spectrum_init.csv", [bin_centers spectrum])

    Nbins = 100

    bins, bin_centers = construct_bins_and_centers(log10(1.e0), log10(1.e3), Nbins)

    bin_centers = 10.0 .^ bin_centers
    bins = 10.0 .^ bins

    dγ = [bins[i+1] - bins[i] for i = 1:length(bins)-1]

    for i = 1:length(shapes), j = 1:length(fieldsolvers)

        sim = joinpath(sim_path, shapes[i] * "/" * fieldsolvers[j])
        info = EntityData(sim)

        γm1 = read_E(info, time, 1, downstr_xmin, downstr_xmax)
        γm1_hist = get_histograms(bins, γm1)
        spectrum = (bin_centers).^2 .* γm1_hist ./ dγ
        writedlm(sim * "/spectrum_downstream_$(@sprintf("%08i", time)).csv", [bin_centers spectrum])

        γm1 = read_E(info, time, 1, upstr_xmin, upstr_xmax)
        γm1_hist = get_histograms(bins, γm1)
        spectrum = (bin_centers).^2 .* γm1_hist ./ dγ
        writedlm(sim * "/spectrum_upstream_$(@sprintf("%08i", time)).csv", [bin_centers spectrum])

    end

end

function Fig19(sim_path, time, plot_name)


    i_colors = [1, 10]
    linestyles = [:dash, :solid]
    shapes = ["zigzag", "esirkepov_3rd"]
    fieldsolvers = ["_Yee", "_min3"]

    f = Figure(size=(1100, 350))



    ax_d = Axis(f[1, 1], xscale=log10, yscale=log10, xlabel=L"\gamma - 1", ylabel=L"\gamma^2 dN/d\gamma",
                        xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9), title = "Downstream")
    xlims!(ax_d, (1.0, 1.e3) )
    ylims!(ax_d, (5.0e-2, 1.e2) )

    ax_u = Axis(f[1, 2], xscale=log10, yscale=log10, xlabel=L"\gamma - 1", ylabel=L"\gamma^2 dN/d\gamma",
                        xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9), title = "Upstream")
    xlims!(ax_u, (1.0, 1.e3) )
    ylims!(ax_u, (5.0e-2, 1.e2) )

    for i = 1:length(shapes), j = 1:length(fieldsolvers)

        sim = joinpath(sim_path, shapes[i] * "/" * fieldsolvers[j])
        d = readdlm(sim * "/spectrum_downstream_$(@sprintf("%08i", time)).csv")
        bin_centers = d[:, 1]
        spectrum = d[:, 2]
        lines!(ax_d, bin_centers, spectrum, 
            linestyle=linestyles[j], linewidth=2.5, 
            color=i_colors[i], colormap = :batlow, colorrange=(0,12))

        d = readdlm(sim * "/spectrum_upstream_$(@sprintf("%08i", time)).csv")
        bin_centers = d[:, 1]
        spectrum = d[:, 2]
        lines!(ax_u, bin_centers, spectrum, 
            linestyle=linestyles[j], linewidth=2.5, 
            color=i_colors[i], colormap = :batlow, colorrange=(0,12))

    end

    d = readdlm(sim_path * "/spectrum_init.csv")
    bin_centers = d[:, 1]
    spectrum = d[:, 2]
    lines!(ax_u, bin_centers, spectrum, 
            linestyle=:dot, linewidth=1.5, 
            color=:black, label=L"f_0")
    
    lines!(ax_u, [0.0], [0.0], linewidth=2.5, color=i_colors[1], colormap = :batlow, colorrange=(0,12), label = "Zigzag")
    lines!(ax_u, [0.0], [0.0], linewidth=2.5, color=i_colors[2], colormap = :batlow, colorrange=(0,12), label = L"\mathrm{Esirkepov} \:\: \mathcal{S}_3")
    lines!(ax_u, [0.0], [0.0], color=:black, linestyle=:dash, linewidth=2.5, label="Yee")
    lines!(ax_u, [0.0], [0.0], color=:black, linestyle=:solid, linewidth=2.5, label="min3")
    axislegend(ax_u, framevisible=false)

    save(plot_name, f)

end

# define plot style
set_light_theme!(fontsize=20)

sim_path = data_path * "shock/"
time = 5000.0
#write_spectra(sim_path, time)

plot_name = plot_path * "Fig19.pdf"
Fig19(sim_path, time, plot_name)