include(joinpath(@__DIR__, "config.jl"))

using CairoMakie, MakieHelper
using EntityTools
using Distributions
using StatsBase
using Printf
using ProgressMeter
using DelimitedFiles


function write_spectra(sim_paths, times)

    cells = ["0-5", "1", 
            "2", "4", "8", "16"] .* "cell"

    titles = [L"\text{Zigzag}", 
            L"\text{Zigzag} + 8 \mathcal{F}", 
            L"\text{Esirkepov} \:\: \mathcal{S}_2", 
            L"\text{Esirkepov} \:\: \mathcal{S}_2 + 8 \mathcal{F}", 
            L"\text{Esirkepov} \:\: \mathcal{S}_3", 
            L"\text{Esirkepov} \:\: \mathcal{S}_5",
            L"\text{Esirkepov} \:\: \mathcal{S}_{9}"]

    labels = [L"\Delta x = 2 d_e",
              L"\Delta x = 1 d_e", 
              L"\Delta x = \frac{1}{2} d_e", 
              L"\Delta x = \frac{1}{4} d_e", 
              L"\Delta x = \frac{1}{8} d_e",
              L"\Delta x = \frac{1}{16} d_e"]

    for i = 1:length(sim_paths), k = 1:length(cells)
        println("reading data from $(sim_paths[i])")
        info = EntityData(sim_paths[i] * "_" * cells[k])

        spectra_1 = read_spectrum(info, 1, t=0.0)
        Nbins = length(spectra_1["E"])-1

            # construct bin centers
        data = Matrix{Float64}(undef, Nbins, length(times)+1)

        for i = 1:Nbins
            data[i, 1] = 10.0^(0.5 * (log10(spectra_1["E"][i]) + log10(spectra_1["E"][i+1])))
        end

        @showprogress for j = 1:length(times)
            spectra_1 = read_spectrum(info, 1, t=times[j], verbose=false)

            for k = 1:Nbins
                data[k, j+1] = spectra_1["N"][k]
            end
        end

        writedlm(sim_paths[i] * "_" * cells[k] * "/spectra.txt", data)
    end

end


function Fig06(sim_paths, times, plot_name)

    t_lim = (0.0, 10000.0)

    #f = Figure(size = (1650, 1000))
    f = Figure(size = (1300, 800))
    
    cmap = :vik
    # colorbars 
    Colorbar(f[1, 2:length(sim_paths)-1], colormap=cmap,  label=L"\mathrm{Time} \:\: t \quad [\omega_p^{-1}]", limits=t_lim, vertical=false, tickalign=1, ticks = 0:1000:10000)

    cells = ["0-5", "1", 
            "2", "4", "8", "16"] .* "cell"

    titles = [L"\text{Zigzag}", 
            L"\text{Zigzag} + 8 \mathcal{F}", 
            L"\text{Esirkepov} \:\: \mathcal{S}_2", 
            L"\text{Esirkepov} \:\: \mathcal{S}_2 + 8 \mathcal{F}", 
            L"\text{Esirkepov} \:\: \mathcal{S}_3", 
            L"\text{Esirkepov} \:\: \mathcal{S}_5",
            L"\text{Esirkepov} \:\: \mathcal{S}_{9}"]

    labels = [L"d_e = \frac{1}{2} \: \Delta x ",
            L"d_e = \Delta x ", 
            L"d_e = 2 \: \Delta x", 
            L"d_e = 4 \: \Delta x", 
            L"d_e = 8 \: \Delta x",
            L"d_e = 16 \: \Delta x"]

    for i = 1:length(sim_paths)

        for k = 1:length(cells)
        println("reading data from $(sim_paths[i])")
        #info = EntityData(sim_paths[i] * "_" * cells[k])
        
        if k == 1
            titlevisible = true
        else
            titlevisible = false 
        end

        ax1 = Axis(f[k+1,i], xlabel=L"\gamma - 1", ylabel=L"\gamma \frac{dN}{d\gamma}",
                xscale=log10, yscale=log10, 
                title=titles[i],
                titlevisible=titlevisible,
                xminorticks = IntervalsBetween(9), yminorticks = IntervalsBetween(9),
                aspect=1.5)

        xlims!(ax1, (1.e-3, 10.0))
        ylims!(ax1, (1.0, 1.e5))

        #spectra_1 = read_spectrum(info, 1, t=0.0)
        #Nbins = length(spectra_1["E"])-1

            # construct bin centers
        # γ = Vector{Float64}(undef, Nbins)

        # for i = 1:Nbins
        #     γ[i] = 0.5 * (spectra_1["E"][i] + spectra_1["E"][i+1])
        # end

        spec = readdlm(sim_paths[i] * "_" * cells[k] * "/spectra.txt")


        @showprogress for j = 1:length(times)
            #spectra_1 = read_spectrum(info, 1, t=times[j], verbose=false)
            
            lines!(ax1, spec[:,1], spec[:,1] .* spec[:,j+1], color = times[j], colormap=cmap, colorrange=t_lim)
        end
        text!(ax1, 1.5e-3, 5e4, text=labels[k], align = (:left, :top), fontsize=15)

        if k < length(cells)
            hidexdecorations!(ax1, ticks = false, minorticks=false)
        end
        
        if i > 1
            hideydecorations!(ax1, ticks = false, minorticks=false)
        end

    end
    end
    resize_to_layout!(f)    
    save(plot_name, f)

end

set_light_theme!()

sim_paths = data_path * "numerical_heating/2D/" .*
            ["zigzag", "zigzag_8filter", 
             "esirkepov_2nd", "esirkepov_2nd_8filter",
             "esirkepov_3rd", 
            "esirkepov_5th", "esirkepov_9th"
                ] .* "/heating"

times = 0:100:10000

#write_spectra(sim_paths, times)

plot_name = plot_path * "Fig06.pdf"
Fig06(sim_paths, times, plot_name)