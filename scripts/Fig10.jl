@info "loading packages"
using CairoMakie,MakieHelper
using MakieExtra
using ColorSchemes, Colors
@info "makie done"
using Printf
using EntityTools
using DelimitedFiles
@info "entitytools done"

function Fig10(sim_paths, run_labels, t, plot_name)

    Nruns = length(sim_paths)
    f = Figure(size=(2100, 600))

    x_l = (-128.0, 128.)
    y_l = x_l
    c_map = Reverse(:RdBu)
    #c_map = :bwr
    myscale = Symlog10(-0.01, 0.01, linscale=1.0)

    jz_lim = (0.01, 1.0)
    jz_lim = (-0.3, 0.3)
    #jz_lim = (-1, 1)

    zoom_center1 = [-100.0, 64.0] ./ 2
    zoom_center2 = [64.0, 96.0] ./ 2
    zoom_width = 32.0 .* [-1.0, 1.0] ./ 2

    Colorbar(f[2, Nruns+1], limits=jz_lim, 
             colormap=c_map, #scale=myscale,
             #ticks = [-1, -0.1, -0.01, 0, 0.01, 0.1, 1],
             ticks = [-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3],
             #minorticksvisible=false,
             label=L"J_z", 
             width=20, height=280
             )

    for i = 1:Nruns

        if i == Nruns 
            x_ticks = [-128, -64, 0, 64, 128]
        else
            x_ticks = [-128, -64, 0, 64]
        end
        @info "reading data"
        info = EntityData(sim_paths[i])
        J = read_field(info, ["J3"], verbose=true; t)
        println(maximum(J["J3"]))
        # Bz
        ax = Axis(f[2, i], 
                ylabel=L"y \:\: \left[ c/\omega_p \right]", 
                xlabel=L"x \:\: \left[ c/\omega_p \right]", 
                aspect = DataAspect(),
                xticks = x_ticks,
                yticks = [-64, 0, 64, 128]
                )
        xlims!(ax, x_l )
        ylims!(ax, y_l )

        image!(ax, x_l, y_l, J["J3"], 
                colormap=c_map, colorrange=jz_lim, #colorscale=myscale,
                rasterize=2)
        if i > 1
            hideydecorations!(ax, ticks = false)
        end

        #text!(ax, (-128+16, -128+16), text = run_labels[i], align = (:left, :center),
        #        offset = (0, 0), color = :black)

        textlabel!(ax, Point2f(-128+4, -128+16), text = run_labels[i], text_align = (:left, :center),
                    offset = (0, 0), text_color = :black, alpha=0.75, strokewidth=0)
        # text!(ax, (-512+64, -512+64), text = run_labels[i], align = (:left, :center),
        #         offset = (0, 0), color = :white)

        ga = f[1, i] = GridLayout()
        colgap!(ga, 0)
        rowgap!(ga, 0)
        ax_zoom1 = Axis(ga[1, 1], 
                aspect = DataAspect())
        xlims!(ax_zoom1, zoom_center1[1] .+ zoom_width )
        ylims!(ax_zoom1, zoom_center1[2] .+ zoom_width )
        heatmap!(ax_zoom1, J["X1"], J["X2"], J["J3"], 
                colormap=c_map, colorrange=jz_lim, #colorscale=myscale,
                rasterize=2)
        hidexdecorations!(ax_zoom1, ticks = false)
        hideydecorations!(ax_zoom1, ticks = false)

        ax_zoom2 = Axis(ga[1, 2], aspect = DataAspect())
        xlims!(ax_zoom2, zoom_center2[1] .+ zoom_width )
        ylims!(ax_zoom2, zoom_center2[2] .+ zoom_width )
        heatmap!(ax_zoom2, J["X1"], J["X2"], J["J3"], 
                colormap=c_map, colorrange=jz_lim, #colorscale=myscale,
                rasterize=2)
        hidexdecorations!(ax_zoom2, ticks = false)
        hideydecorations!(ax_zoom2, ticks = false)

        rowsize!(f.layout, 1, Aspect(1, 0.5))
        zoom_lines!(ax, ax_zoom1,
                rectattrs=(linestyle=:solid, strokecolor=:white, strokewidth=3.0, transparency=0.5),
                lineattrs=(linestyle=:solid, color=:white, linewidth=3.0, transparency=0.5))
        zoom_lines!(ax, ax_zoom1,
                rectattrs=(linestyle=:dash, strokewidth=2.0, transparency=0.5),
                lineattrs=(linestyle=:dash, linewidth=2.0, transparency=0.5))
                
        zoom_lines!(ax, ax_zoom2,
                rectattrs=(linestyle=:solid, strokecolor=:white, strokewidth=3.0, transparency=0.5),
                lineattrs=(linestyle=:solid, color=:white, linewidth=3.0, transparency=0.5))
        zoom_lines!(ax, ax_zoom2,
                rectattrs=(linestyle=:dash, strokewidth=2.0, transparency=0.5),
                lineattrs=(linestyle=:dash, linewidth=2.0, transparency=0.5))

    end

    resize_to_layout!(f)
    @info "saving"
    save(plot_name, f)

    @info "done"
    GC.gc()
end

set_light_theme!(fontsize=30)

runs = ["Zigzag", 
        L"\text{Esirkepov} \:\: \mathcal{S}_2", 
        L"\text{Esirkepov} \:\: \mathcal{S}_3", 
        L"\text{Esirkepov} \:\: \mathcal{S}_5",
        L"\text{Esirkepov} \:\: \mathcal{S}_9",
        L"\text{Zigzag} \: + 16 \: \mathcal{F}"]

sim_paths = data_path * "decaying_turb/2D/" .*
            ["zigzag", "esirkepov_2nd", "esirkepov_3rd", "esirkepov_5th", "esirkepov_9th", "zigzag_16filters"] .* "/turbulence"

t = 1092.0

plot_name = plot_path * "Fig10.pdf"
Fig10(sim_paths, runs, t, plot_name)
