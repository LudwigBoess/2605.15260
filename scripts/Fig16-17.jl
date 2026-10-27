include(joinpath(@__DIR__, "config.jl"))

using CairoMakie, MakieHelper
using EntityTools
using StatsBase
using Printf
using ProgressMeter
using DelimitedFiles

abs_sum(x) = sum(abs.(x))
abs_mean(x) = sum(abs.(x)) / length(x)

function Fig_16(sim_path, deposits, t, plot_name)

    deposit_names = [L"\text{Zigzag} \: + \: 5 \mathcal{F}", L"\mathrm{Esirkepov} \: \mathcal{S}_3 \: + \: 5 \mathcal{F}"]
    runs = ["Yee", "Cowan", "Lehe", "min1", "min2", "min3", 
            "min4", "min5", "min6"]
    arr_labels = [L"\textbf{Yee}", L"\textbf{Cowan}", L"\textbf{Lehe}\leftrightarrow", 
            L"\textbf{min1}", L"\textbf{min2}", L"\textbf{min3}", L"\textbf{min4}", L"\textbf{min5}", L"\textbf{min6}\leftrightarrow"]
    Nruns = length(runs)
    println("Nruns: ", Nruns)

    cmap = :seaborn_colorblind
    Nmax = Nruns

    f = Figure(size=(850, 900))
    
    for j = 1:2
        ax_Bz = Axis(f[1,j], xlabel = L"t \quad [10^3 \: \omega_p^{-1}]", ylabel = L"\langle B_z \rangle \:\: \left[ \sqrt{4\pi n_e m_e c^2} \right]",
                    yscale = log10, yminorticks = IntervalsBetween(9), title = deposit_names[j], aspect=1,
                    xticks = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
                    )
        xlims!(ax_Bz, 0.0, 10.0)
        ylims!(ax_Bz, 1.e-5, 1.0)
        hidexdecorations!(ax_Bz, ticks = false, minorticks=false)
        if j == 2
            hideydecorations!(ax_Bz, ticks = false, minorticks=false)
        end
        ax_T = Axis(f[2,j], xlabel = L"t \quad [10^3 \: \omega_p^{-1}]", ylabel = L"\log_{10} \: (T - T_0)/T_0",
                    yticks = -6:3, aspect=1, 
                    xticks = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
                    )
        xlims!(ax_T, 0.0, 10.0)
        ylims!(ax_T, -6, 3)
        if j == 2
            hideydecorations!(ax_T, ticks = false, minorticks=false)
        end
        for i = 1:Nruns
            @info "reading data: $(runs[i]) $(deposits[j])"
            d = readdlm(sim_path * "/data/heating_new2_$(deposits[j])_$(runs[i]).txt")
            data = Dict("t" => d[:,1], "B3" => d[:,2], "T22" => abs.(d[:,3] .- d[1,3]) ./ d[1,3])  # normalize temperature increase
            lines!(ax_Bz, 1.e-3 .* data["t"], data["B3"], color = i , colormap = cmap, colorrange=(1,Nmax), 
                        linewidth=4, label=runs[i])
            lines!(ax_T, 1.e-3 .* data["t"], log10.(data["T22"]), color = i , colormap = cmap, colorrange=(1,Nmax), 
                    linewidth=4)
        end
        if j == 2
            lins = [lines!(ax_Bz, [0.0], [0.0], color = i, colormap = cmap, colorrange=(1,Nmax), linewidth=3) for i = 1:Nruns]
            Legend(f[3, 1:2], lins, arr_labels, framevisible=false, orientation = :horizontal, nbanks=2)
        end
    end

    @info "saving"
    save(plot_name, f)

    @info "done"
    # de-allocate data (just for caution)
    GC.gc()
end

set_light_theme!(fontsize=25)

sim_path = data_path * "drifting_plasma/2D/"
deposits = ["fieldstencil_zigzag", "fieldstencil_e3"]
t = 10000.0
plot_name = plot_path * "Fig16.pdf"
Fig_16(sim_path, deposits, t, plot_name)


function Fig_17(sim_path, deposits, t, plot_name)

    deposit_names = [L"\text{Zigzag} \: + \: 5 \mathcal{F}", L"\mathrm{Esirkepov} \: \mathcal{S}_3 \: + \: 5 \mathcal{F}"]
    runs = ["Yee", "Lehe", "min3D", "min1D-1", "min1D-2"]
    arr_labels = [L"\textbf{Yee}", L"\textbf{Lehe3D}\leftrightarrow", 
            L"\textbf{min3D1}", L"\textbf{min3D2}", L"\textbf{min3D3}"]
    Nruns = length(runs)

    cmap = :seaborn_colorblind
    Nmax = 9

    color_idxs = [1, 2, 3, 4, 6, 7, 8, 9]

    f = Figure(size=(870, 900))
    
    for j = 1:2
        ax_Bz = Axis(f[1,j], xlabel = L"t \quad [10^3 \: \omega_p^{-1}]", ylabel = L"\langle B_z \rangle \:\: \left[ \sqrt{4\pi n_e m_e c^2} \right]",
                    yscale = log10, yminorticks = IntervalsBetween(9), title = deposit_names[j], aspect=1,
                    xticks = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
                    )
        xlims!(ax_Bz, 0.0, 10.0)
        ylims!(ax_Bz, 1.e-7, 1.0)
        hidexdecorations!(ax_Bz, ticks = false, minorticks=false)
        if j == 2
            hideydecorations!(ax_Bz, ticks = false, minorticks=false)
        end
        ax_T = Axis(f[2,j], xlabel = L"t \quad [10^3 \: \omega_p^{-1}]", ylabel = L"\log_{10} \: (T - T_0)/T_0",
                    yticks = -6:3, aspect=1, 
                    xticks = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
                    )
        xlims!(ax_T, 0.0, 10.0)
        ylims!(ax_T, -6, 3)
        if j == 2
            hideydecorations!(ax_T, ticks = false, minorticks=false)
        end
        for i = 1:Nruns
            file_path = sim_path * deposits[j] * "/_$(runs[i])/_$(runs[i])_stats.csv"
            if !isfile(file_path)
                @warn "file $file_path does not exist, skipping"
                continue
            end
            @info "reading data: $(runs[i]) $(deposits[j])"
            d = readdlm(file_path, ',', skipstart=1)
            data = Dict("t" => d[:,2], "B3" => d[:,5], "T22" => abs.(d[:,6] .- d[1,6]) ./ d[1,6])  # normalize temperature increase

            lines!(ax_Bz, 1.e-3 .* data["t"], data["B3"], color = color_idxs[i], colormap = cmap, colorrange=(1,Nmax), 
                        linewidth=4, label=runs[i])
            lines!(ax_T, 1.e-3 .* data["t"], log10.(data["T22"]), color = color_idxs[i], colormap = cmap, colorrange=(1,Nmax), 
                    linewidth=4)
        end
        if j == 2
            lins = [lines!(ax_Bz, [0.0], [0.0], color = color_idxs[i], colormap = cmap, colorrange=(1,Nmax), linewidth=3) for i = 1:Nruns]
            Legend(f[3, 1:2], lins, arr_labels, framevisible=false, orientation = :horizontal, nbanks=1)
        end
    end

    #axislegend(ax_Bz, framevisible=false)
    @info "saving"
    save(plot_name, f)

    @info "done"
    # de-allocate data (just for caution)
    GC.gc()
end


set_light_theme!(fontsize=25)

sim_path = data_path * "drifting_plasma/3D/"
deposits = ["zigzag", "esirkepov_3rd"]
t = 10000.0
plot_name = plot_path * "Fig17.pdf"
Fig_17(sim_path, deposits, t, plot_name)
