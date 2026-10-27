include(joinpath(@__DIR__, "config.jl"))

using Printf
using CairoMakie, MakieHelper
using EntityTools

function read_Blinne_data(fi)
    B3 = h5read(fi, "/Step0/fB3")

    X = h5read(fi, "/Step0/X1e") .* 10.0
    Y = h5read(fi, "/Step0/X2e") .* 10.0

    return X, Y, B3
end

function read_particle_data(fi)
    X = h5read(fi, "/Step0/pX1_1")[1] .* 10.0
    Y = h5read(fi, "/Step0/pX2_1")[1] .* 10.0
    return X, Y
end


function Fig14(sim_path, arr_labels, run_names, plot_name)

    Nruns = length(run_names)
    f = Figure(size=(2000, 700))

    x_l = (-0.1, 0.1)
    y_l = (-0.1, 0.1)
    offset = 0.1
    b_lim = (-10.0, 10.0)
    c_map = :RdBu

    deg_labels = ["0", "45", "90"]
    myscale = Symlog10(0.05, linscale=0.5)

    Colorbar(f[1:3, Nruns+1], limits=b_lim, colormap=c_map, label=L"B_z \:\: [a.u.]", vertical=true, tickalign=1,
        scale=myscale, ticks=[-10, -1, -0.1, -0.05, 0, 0.05,  0.1, 1, 10], minorticksvisible=false)

    for j = 1:3, i = 1:Nruns
        @info "plotting $(run_names[i])_$(deg_labels[j])deg"
        out_path = sim_path * "_$(run_names[i])_$(deg_labels[j])deg"
        info = EntityData(out_path)

        data = read_field(info, ["B3"], t=0.45)
        p = read_particles(info, 1, ["X1", "X2"], t=0.45)
        
        title_is_visible = false
        if j == 1
            title_is_visible = true
        end

        ax = Axis(f[j, i], ylabel=L"\Delta y", xlabel=L"\Delta x", aspect = DataAspect(),
                title = arr_labels[i], titlevisible=title_is_visible, xticks=[-0.05, 0.05], yticks=[-0.05, 0.05])
        xlims!(ax, x_l)
        ylims!(ax, y_l)

        heatmap!(ax, data["X1"] .- p["X1"], data["X2"] .- p["X2"], data["B3"] .* 1.e6, colormap=c_map, colorrange=b_lim,
            colorscale=myscale, rasterize=2)

        colgap!(f.layout, i, 15)
        if j < 3
            rowgap!(f.layout, j, 15)
        end

        if j < 3
            hidexdecorations!(ax, ticks = false)
        end        
        if i > 1
            hideydecorations!(ax, ticks = false)
        end

    end

    @info "saving"
    save(plot_name, f)

    @info "done"
    GC.gc()
end

set_light_theme!(fontsize=30)
sim_path = data_path * "single_particle/"
run_names = ["Yee", "Cowan", "Lehe", "min1", "min2", "min3", "min4", "min5", "min6"]
arr_labels = [L"\textbf{Yee}", L"\textbf{Cowan}", L"\textbf{Lehe}\leftrightarrow", 
            L"\textbf{min1}", L"\textbf{min2}", L"\textbf{min3}", L"\textbf{min4}", L"\textbf{min5}", L"\textbf{min6}\leftrightarrow"]
plot_name = plot_path * "Fig14.pdf"

Fig14(sim_path, arr_labels, run_names, plot_name)