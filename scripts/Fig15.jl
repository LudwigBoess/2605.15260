include(joinpath(@__DIR__, "config.jl"))

using CairoMakie, MakieHelper
using Printf
using StatsBase
using EntityTools
using DelimitedFiles


function Fig15(sim_path, deposits, t, plot_name)

    runs = ["Yee", "Cowan", "Lehe", "min1", "min2", "min3", "min4", "min5", "min6"]
    arr_labels = [L"\textbf{Yee}", L"\textbf{Cowan}", L"\textbf{Lehe}\leftrightarrow", 
            L"\textbf{min1}", L"\textbf{min2}", L"\textbf{min3}", L"\textbf{min4}", L"\textbf{min5}", L"\textbf{min6}\leftrightarrow"]
    Nruns = length(runs)

    f = Figure(size=(1150, 830))
    Bzgrid = f[1, 1] = GridLayout()
    Bfiltergrid = f[2, 1] = GridLayout()
    Bkgrid = f[3, 1] = GridLayout()
    N_lim = (0.5, 1.5)
    Bz_lim = (-1.0, 1.0)
    Bk_lim = (1.e-4, 1.e2)
    myscale = Symlog10(1.e-4, linscale=1.e-4)
    myticks=[-1, -0.1, -0.01, -0.001, 0.0, 0.001, 0.01, 0.1, 1]
    deposit_names = [L"\text{Zigzag}", L"\text{Esirkepov} \: \mathcal{S}_3"]

    for i = 1:Nruns

        for j = 1:2
            @info "reading data: $(runs[i]) $(deposits[j])"

            if j == 1
                titlevisible = true
            else
                titlevisible = false 
            end
            # read one file to check if it exists
            info = EntityData(sim_path * deposits[j] * "/_" * runs[i])
            data = read_field(info, ["N", "B3"], verbose=true; t)

            # Bz
            @info "plotting Bz"
            ax = Axis(Bzgrid[j,i], ylabel=L"y \:\: \left[ c/\omega_e \right]", xlabel=L"x \:\: \left[ c/\omega_e \right]", 
            xticks = LinearTicks(3), yticks = LinearTicks(3), aspect = DataAspect(),
            title = arr_labels[i], titlevisible=titlevisible)
            heatmap!(ax, data["X1"], data["X2"], data["B3"],
                    colormap=:RdBu, colorrange=Bz_lim, colorscale=myscale,
                    rasterize=2)
            if i > 1
                hideydecorations!(ax, ticks = false)
            end
            if j == 1
                hidexdecorations!(ax, ticks = false)
            end
            if i == 1
                text!(ax, -4.5, 4.5, text=deposit_names[j], align = (:left, :top), fontsize=15, color=:white)
            end
            

            @info "plotting k-space"
            # N cells per skin depth
            Ncells = 256/10

            k = fftshift(fftfreq(size(data["B3"], 1), 1/Ncells))
            k_B3 = fftshift(fft(data["B3"]))
            ftt_x_l = 100 .* (minimum(k), maximum(k))
            ftt_y_l = 100 .* (minimum(k), maximum(k))

            # FFT
            ax = Axis(Bkgrid[j,i], ylabel=L"k_y \:\: \left[ 10^2 \: \omega_e/c \right]", xlabel=L"k_x \:\: \left[ 10^2 \: \omega_e/c \right]", 
            xticks = LinearTicks(3), yticks = LinearTicks(3), aspect = DataAspect())

            fft_image = abs.(real(k_B3))
            image!(ax, ftt_x_l, ftt_y_l, fft_image, colormap=:haline, colorrange=Bk_lim, colorscale=log10,
                    rasterize=2)
            
            if i > 1
                hideydecorations!(ax, ticks = false)
            end
            if j == 1
                hidexdecorations!(ax, ticks = false)
            end
            if i == 1
                text!(ax, -1.8, 1.8, text=deposit_names[j], align = (:left, :top), fontsize=15, color=:white)
            end

            # filtered image
            @info "plotting filtered Bz"
            k_B3 = fft(data["B3"])

            ax = Axis(Bfiltergrid[j,i], ylabel=L"y \:\: \left[ c/\omega_e \right]", xlabel=L"x \:\: \left[ c/\omega_e \right]", 
            xticks = LinearTicks(3), yticks = LinearTicks(3), aspect = DataAspect())

            fft_image = abs.(real(k_B3))
            sel = findall( fft_image .< 0.9 * maximum(fft_image) )
            k_B3[sel] .= 0.0
            println(size(k_B3))

            B_image = real(ifft(k_B3))
            println(size(B_image))
            println(size(B_image))

            heatmap!(ax, data["X1"], data["X2"], B_image,
                    colormap=:RdBu, colorrange=Bz_lim, colorscale=myscale,
                    rasterize=2)
            
            if i > 1
                hideydecorations!(ax, ticks = false)
            end
            if j == 1
                hidexdecorations!(ax, ticks = false)
            end
            if i == 1
                text!(ax, -4.5, 4.5, text=deposit_names[j], align = (:left, :top), fontsize=15, color=:white)
            end
        end

    end

    Colorbar(f[1, 2], limits=Bz_lim, colormap=:RdBu, label=L"B_z(x) \:\: \left[ \sqrt{4\pi n_e m_e c^2} \right]", tickalign=1, scale=myscale, ticks = myticks, minorticks=IntervalsBetween(9),
    tickformat = values -> [iszero(value) ? "0" : rich("$(value < 0 ? "-" : "")10", superscript("$(@sprintf("%i", log10(abs(value))))")) for value in values])
    colgap!(Bzgrid, 15.0)
    rowgap!(Bzgrid, 15.0)

    Colorbar(f[2, 2], limits=Bz_lim, colormap=:RdBu, label=L"B_{z, k_{90\%}}(x) \:\: \left[ \sqrt{4\pi n_e m_e c^2} \right]", tickalign=1, scale=myscale, ticks = myticks, minorticks=IntervalsBetween(9),
    tickformat = values -> [iszero(value) ? "0" : rich("$(value < 0 ? "-" : "")10", superscript("$(@sprintf("%i", log10(abs(value))))")) for value in values])
    colgap!(Bfiltergrid, 15.0)
    rowgap!(Bfiltergrid, 15.0)

    Colorbar(f[3, 2], limits=Bk_lim, colormap=:haline, label=L"|B_z(\hat{k})| \:\: \left[ a.u. \right]", tickalign=1, scale=log10, minorticks=IntervalsBetween(9),
        tickformat = values -> [rich("$(value < 0 ? "-" : "")10", superscript("$(@sprintf("%i", log10(abs(value))))")) for value in values])

    colgap!(Bkgrid, 15.0)
    rowgap!(Bkgrid, 15.0)

    @info "saving"
    save(plot_name, f)

    @info "done"
    # de-allocate data (just for caution)
    GC.gc()
end

set_light_theme!(fontsize=18)

sim_path = data_path * "drifting_plasma/2D/"
deposits = ["fieldstencil_zigzag", "fieldstencil_e3"]
t = 10000.0
plot_name = plot_path * "Fig15.pdf"
Fig15(sim_path, deposits, t, plot_name)
