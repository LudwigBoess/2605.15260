include(joinpath(@__DIR__, "config.jl"))

using CairoMakie, MakieHelper
using Printf
using StatsBase
using EntityTools
using DelimitedFiles


function Fig18(sim_path, deposits, t, plot_name)

    f = Figure(size=(1250, 920))
    Ngrid = f[1, 1] = GridLayout()
    Tgrid = f[2, 1] = GridLayout()
    Bzgrid = f[3, 1] = GridLayout()

    N_lim = (0.0, 4.0)
    Bz_lim = (-2.0, 2.0)
    T_lim = (1.0e-4, 10.0)

    x_l = (1_500.0, 4_500.0)
    y_l = (-387.5, 387.5)

    runs = ["Yee", "min3"]
    deposit_names = [L"\text{Zigzag} \: + \: 8 \mathcal{F}", L"\text{Esirkepov} \: \mathcal{S}_3 + \: 8 \mathcal{F}"]
    Nruns = length(runs)

    for i = 1:Nruns

        for j = 1:2
            @info "reading data: $(runs[i]) $(deposits[j])"

            if i == 1
                titlevisible = true
            else
                titlevisible = false 
            end
            # read one file to check if it exists
            info = EntityData(sim_path * deposits[j] * "/_" * runs[i])
            data = read_field(info, ["N_1", "N_2", "B3", "T22_1"], verbose=true, slice = [ [x_l[1], x_l[2]], [y_l[1], y_l[2]] ]; t)
            
            Temp = data["T22_1"] ./ data["N_1"]
            Temp[Temp .<= 0.0] .= 1.e-4 # set negative temperatures to NaN for better plotting
            Temp[isnan.(Temp)] .= 1.e-4 # set NaN temperatures to a small value for better plotting
            
            println("minmax(Temp) = ", extrema(Temp))

            ax = Axis(Ngrid[i,j], ylabel=L"y \:\: \left[ c/\omega_e \right]", xlabel=L"x \:\: \left[ c/\omega_e \right]", 
                    aspect = DataAspect(),
                    title = deposit_names[j], titlevisible=titlevisible)
            xlims!(ax, x_l )
            ylims!(ax, y_l )
            image!(ax, x_l, y_l, data["N_1"] .+ data["N_2"],
                    colormap=:viridis, colorrange=N_lim,
                    rasterize=2)

            hidexdecorations!(ax, ticks = false, minorticks = false)

            if j == 2
                hideydecorations!(ax, ticks = false)
            end
            text!(ax, 1600.0, 350.0, text=runs[i], align = (:left, :top), fontsize=25, color=:black)

            ax = Axis(Tgrid[i,j], ylabel=L"y \:\: \left[ c/\omega_e \right]", xlabel=L"x \:\: \left[ c/\omega_e \right]", 
                    aspect = DataAspect())
            xlims!(ax, x_l )
            ylims!(ax, y_l )
            image!(ax, x_l, y_l, Temp,
                    colormap=:thermal, colorrange=T_lim,
                    colorscale=log10,
                    rasterize=2)

            hidexdecorations!(ax, ticks = false, minorticks = false)

            if j == 2
                hideydecorations!(ax, ticks = false)
            end
            text!(ax, 1600.0, 350.0, text=runs[i], align = (:left, :top), fontsize=25, color=:black)


            ax = Axis(Bzgrid[i,j], ylabel=L"y \:\: \left[ c/\omega_e \right]", xlabel=L"x \:\: \left[ c/\omega_e \right]", 
                    aspect = DataAspect())
            xlims!(ax, x_l )
            ylims!(ax, y_l )
            image!(ax, x_l, y_l, data["B3"],
                    colormap=:RdBu, colorrange=Bz_lim,
                    rasterize=2)

            if i < Nruns
                hidexdecorations!(ax, ticks = false, minorticks = false)
            end
            if j == 2
                hideydecorations!(ax, ticks = false)
            end
            text!(ax, 1600.0, 350.0, text=runs[i], align = (:left, :top), fontsize=25, color=:black)

        end
    end

    Colorbar(f[1, 2], limits=N_lim, colormap=:viridis, label=L"n_\pm \:\: \left[ n_{0} \right]", tickalign=1)
    colgap!(Ngrid, 15.0)
    rowgap!(Ngrid, 15.0)

    Colorbar(f[2, 2], limits=T_lim, colormap=:thermal, 
                scale = log10, minorticks = IntervalsBetween(9), 
                label=L"T \:\: \left[ m_\pm c^2 \right]", tickalign=1)
    colgap!(Tgrid, 15.0)
    rowgap!(Tgrid, 15.0)

    Colorbar(f[3, 2], limits=Bz_lim, colormap=:RdBu, label=L"B_z \:\: \left[ \sqrt{4\pi n_\pm m_\pm c^2} \right]", tickalign=1)
    colgap!(Bzgrid, 15.0)
    rowgap!(Bzgrid, 15.0)

    resize_to_layout!(f)

    @info "saving"
    save(plot_name, f)

    @info "done"
    # de-allocate data (just for caution)
    GC.gc()

end

set_light_theme!(fontsize=20)

sim_path = data_path * "shock/"
deposits = ["zigzag", "esirkepov_3rd"]
t = 5000.0
plot_name = plot_path * "Fig18.pdf"
Fig18(sim_path, deposits, t, plot_name)
