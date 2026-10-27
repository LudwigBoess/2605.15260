include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using CairoMakie, MakieHelper
using ColorSchemes
@info "makie done"
using Printf
using EntityTools
using DelimitedFiles
@info "entitytools done"

function Fig08(sim_paths, run_labels, t, plot_name)

    Nruns = length(sim_paths)
    f = Figure(size=(1200, 500))

    x_l = (-12.8, 12.8)
    y_l = x_l
    c_map = :RdBu
    jz_lim = (-0.01, 0.01)
    myscale = Symlog10(0.002, linscale=1.0)


    Colorbar(f[1:2, Nruns+1], limits=jz_lim, 
             colormap=c_map, scale=myscale,
             label=L"E_y \:\: \left[ \sqrt{8 \pi n_b m_e c^2} \right]", 
             width=20, height=280,
             ticks =[-0.01, -0.005, -0.002, -0.001, 0, 0.001, 0.002, 0.005, 0.01], minorticksvisible=false
             )

    #x_ticks = [-10.0, -5.0, 0.0, 5.0, 10.0]
    x_ticks = [-8, 0, 8]

    subfolder = ["", "_8filters"]
    subfolder_labels = [L"0 \: \mathcal{F}", L"8 \: \mathcal{F}"]

    for i = 1:Nruns, j = 1:2

        @info "reading data"
        info = EntityData(sim_paths[i] * subfolder[j])
        d = read_field(info, ["E2"], verbose=true; t)
        # Bz
        ax = Axis(f[j, i], 
                ylabel=L"y \:\: \left[ c/\omega_p \right]", 
                xlabel=L"x \:\: \left[ c/\omega_p \right]", 
                aspect = DataAspect(),
                xticks = x_ticks,
                yticks = x_ticks,
                title = run_labels[i],
                titlevisible = j==1 ? true : false
                )
        xlims!(ax, x_l )
        ylims!(ax, y_l )

        image!(ax, x_l, y_l, d["E2"], 
                colormap=c_map, colorrange=jz_lim, colorscale=myscale,
                rasterize=2)
        if i > 1
            hideydecorations!(ax, ticks = false)
        else
            #text!(ax, (-10, 10), text = subfolder_labels[j], align = (:left, :center),
            #        offset = (0, 0), color = :black)
            textlabel!(ax, Point2f(-8, 10), text = subfolder_labels[j],
                    offset = (0, 0), text_color = :black, alpha=0.75, strokewidth=0)
        end
        if j == 1
            hidexdecorations!(ax, ticks = false)
        end
    end
    colgap!(f.layout, 5)
    rowgap!(f.layout, 5)
    colgap!(f.layout, Nruns, 10)
    resize_to_layout!(f)
    @info "saving"
    save(plot_name, f)

    @info "done"
    GC.gc()
end

set_light_theme!(fontsize=28)

runs = ["Zigzag", 
        L"\text{Esirkepov} \:\: \mathcal{S}_2", 
        L"\text{Esirkepov} \:\: \mathcal{S}_3",
        L"\text{Esirkepov} \:\: \mathcal{S}_5",
        L"\text{Esirkepov} \:\: \mathcal{S}_9"
        ]

sim_paths = data_path * "oblique_modes/" .*
            ["zigzag", 
            "esirkepov_2nd", 
            "esirkepov_3rd", 
            "esirkepov_5th", 
            "esirkepov_9th",
            ] .* "/oblique_modes"

t = 1500.0

plot_name = plot_path * "Fig08.pdf"
Fig08(sim_paths, runs, t, plot_name)
