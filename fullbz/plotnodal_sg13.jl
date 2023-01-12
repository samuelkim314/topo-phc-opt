using Crystalline, DelimitedFiles
using Brillouin, Interpolations
using MPBUtils
include("bz_utils.jl")

makie = false
if makie
    using GLMakie, CairoMakie
    # GLMakie.activate!
else
    using PlotlyJS
end

sgnum    = 13
subdir = "sg$sgnum-dim$D-n4.00/opt-local-4-579-best"
calcname = "$subdir/3944-bzsymreduced"
# subdir = "sg$sgnum-dim$D-n4.00/opt-local-7-1242-best"
# calcname = "$subdir/3918-bzsymreduced"
# subdir = "sg$sgnum-dim$D-n4.00/opt-local-9-1916-best"
# calcname = "$subdir/4417-bzsymreduced"
outname  = "output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
Rs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("input/$calcname-temp.sh") # Rs, flat, isoval, epsin, epsout, kvecs
Gs = reciprocalbasis(Rs)
#plot(flat, Rs; isoval=isoval)

kp = irrfbz_path(sgnum, Rs)

freqs = sort!(data[:,6:end], dims=2)
kvs = collect(eachrow(data[:,2:4]))
kxs, kys, kzs = sort(unique(getindex.(kvs, 1))), sort(unique(getindex.(kvs, 2))), sort(unique(getindex.(kvs, 3)))
Δx, Δy, Δz = -(-)(extrema(kxs)...), -(-)(extrema(kys)...), -(-)(extrema(kzs)...)
Nkx, Nky, Nkz = length(kxs), length(kys), length(kzs)
freqs = reshape(freqs, (Nkx, Nky, Nkz, size(freqs, 2)))

# Extract centroid of the minimum contour in each kz slice
coords = zeros(Nky, 3)
for i in eachindex(kys)
    freqxz = abs.(freqs[:, i, :, 3] - freqs[:, i, :, 2])

    ci = calc2dcentroid(freqxz, kxs, kzs, threshold=2.0)
    cx, cz = ci[1], ci[2]

    coords[i, :] = [cx, kys[i], cz]
end

# Get frequency of the nodal line
freqnodal = zeros(length(kys))
for i in eachindex(kys)
    itp2 = LinearInterpolation((kxs, kys), freqs[:, i, :, 2])
    itp3 = LinearInterpolation((kxs, kys), freqs[:, i, :, 3])
    kxi, kzi = coords[i, 1], coords[i, 3]
    freqiavg = (itp2[kxi, kzi] + itp3[kxi, kzi]) / 2
    freqnodal[i] = freqiavg
end
# p = PlotlyJS.plot(PlotlyJS.scatter(;x=kys, y=freqnodal, mode="line"),
#     Layout(;xaxis_title="k_y", yaxis_title="Frequency", xaxis_showgrid=false))

coords = coords .+ reshape([-1 0 0], 1, 3)  # Moving it to be inside Wigner-Seitz cell 4-579-3944
# coords = coords .+ reshape([0 0 -1], 1, 3)  # Moving it to be inside Wigner-Seitz cell
# for (i, coord) in enumerate(eachrow(coords))
#     if coord[1] < 0
#         coords[i, :] = coord .* [-1, 1, -1] .+ [0, 0, 1]
#         # println((i, coord))
#     else
#         # println((i, coord))
#     end
# end
# coords = coords .+ reshape([0 0 -1], 1, 3)  # Moving it to be inside Wigner-Seitz cell

# Take out some of the coordinates 7-1242-3918
# coords = coords[vcat(1:30, 34:end), :]
# freqnodal = freqnodal[vcat(1:30, 34:end)]

# # Apply symmetry operations to extend to full BZ
coords_list = [coords]
pg = pointgroup(primitivize(spacegroup(sgnum, Val(3))))
for op in pg
    global coords
    rotation(op) == I && continue
    push!(coords_list, (rotation(op)'\coords')')
end
# append!(coords_list, [c .+ reshape([1 0 0], 1, 3) for c in coords_list])
# append!(coords_list, [c .+ reshape([-1 0 0], 1, 3) for c in coords_list])
# append!(coords_list, [c .+ reshape([0 0 1], 1, 3) for c in coords_list])
# append!(coords_list, [c .+ reshape([0 0 -1], 1, 3) for c in coords_list])

# Convert nodal line to Cartesian basis
coords_cart = []
for coordsi in coords_list
    coordsi = [Array(Brillouin.cartesianize(coordi, Gs)) for coordi in eachrow(coordsi)]
    coordsi = hcat(coordsi...)'
    push!(coords_cart, coordsi)
end
coords_list = coords_cart

if makie
    # Plot BZ in Cartesian basis
    c = wignerseitz(Gs)
    fig, ax, p = GLMakie.plot(c)
    # Plot high-symmetry k-path
    GLMakie.plot!(kp, textkws=(; strokewidth=0))

    cm = cgrad([:green1, :red])

    # Plot nodal lines
    nodalplot = nothing
    for coordsi in coords_list
        global nodalplot
        nodalplot = GLMakie.lines!(ax, coordsi[:, 1], coordsi[:, 2], coordsi[:, 3]; color=freqnodal, linewidth=10, 
            # colormap=:diverging_linear_bjr_30_55_c53_n256,
            # colormap=:brg, colorrange=(0.315462299523985, 0.35))
            # colormap=:brg)
            colormap=cm)
    end
    # limits: (0.315462299523985, 0.32873040920883245)
    GLMakie.Colorbar(fig[1, 2], nodalplot)
    # Colorbar(fig[1, 2], nodalplot, limits=(0.315462299523985, 0.32873040920883245))
    # Colorbar(fig[1, 2], nodalplot, colorrange=(0.315462299523985, 0.35))
    # ax.azimuth = -0.1pi
    # ax.elevation = 0.25pi
    # ax.azimuth = 0.4pi
    # ax.elevation = -0.2pi
    display(fig)
    # CairoMakie.save("figures/sg13-nodal.pdf", fig)
else
    # Plot BZ in Cartesian basis
    c = wignerseitz(Gs)
    p = Brillouin.plot(c)
    BZtrace = p.plot.data
    layout = p.plot.layout
    # println(BZtrace)

    # Get rid of red axes
    BZtrace = filter(trace -> trace["hovertext"] == "Cell", BZtrace)

    # Plot HS k-lines (path from irrfbz_path)
    p_k = Brillouin.plot(kp)
    append!(BZtrace, p_k.plot.data)

    # Plot nodal lines
    # nodaltraces = [PlotlyJS.scatter(x=coordsi[:, 1], y=coordsi[:, 2], z=coordsi[:, 3], 
    #     type="scatter3d", mode="lines", marker_color=coordsi[:, 2], line = attr(width=10)) 
    #     for coordsi in coords_list]
    nodaltraces = [PlotlyJS.scatter(x=coordsi[:, 1], y=coordsi[:, 2], z=coordsi[:, 3], 
        type="scatter3d", mode="lines", 
        line=attr(width=10, colorscale="Bluered", color=freqnodal, showscale=true)) 
        for coordsi in coords_list]

    # PlotlyJS.plot([PlotlyJS.scatter(x=coords[:, 1], y=coords[:, 2], z=coords[:, 3], type="scatter3d", mode="lines", marker_color="black")],
    #     Layout(scene_xaxis_range=[-0.5, 0.5], scene_yaxis_range=[-0.5, 0.5], scene_zaxis_range=[-0.5, 0.5]))
    # PlotlyJS.plot(nodaltraces,
    #     Layout(scene_xaxis_range=[-0.5, 0.5], scene_yaxis_range=[-0.5, 0.5], scene_zaxis_range=[-0.5, 0.5]))

    # PlotlyJS.plot(nodaltraces)
    p = PlotlyJS.plot(append!(nodaltraces, BZtrace),
        p.plot.layout)

    # PlotlyJS.savefig(p, "sg13BZ.svg")
end