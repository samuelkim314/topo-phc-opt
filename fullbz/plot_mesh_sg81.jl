using Crystalline, PyPlot, DelimitedFiles, LinearAlgebra, Interpolations, Brillouin
using MPBUtils
using Main.KPaths

sgnum    = 81
D        = 3
res      = 32
subdir = "sg$sgnum-dim$D-n4.00/gaps-ff-direct/rand-local-10-4822-best"
calcname = "$subdir/4358-bzsymreduced"
nbandsbelow = 7
# subdir = "sg$sgnum-dim$D-n4.00/gaps-ff-direct/rand-local-9-3098-best"
# calcname = "$subdir/4990-bzsymreduced"
# nbandsbelow = 7
subdir = "sg$sgnum-dim$D-n4.00/gaps-ff-direct/local-9-3098-best"
calcname = "$subdir/4384-bzsymreduced"
nbandsbelow = 4
subdir = "sg$sgnum-dim$D-n4.00/gaps-ff-direct/local-7-349-best"
calcname = "$subdir/2840-bzsymreduced"
nbandsbelow = 4
outname  = "output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
pRs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("output/$calcname-temp.sh") # Rs, flat, isoval, epsin, epsout, kvecs
cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)
Rs = conventionalize(pRs, cntr)
Gs = reciprocalbasis(pRs)
#plot(flat, Rs; isoval=isoval)
println((Rs, Gs))

# extract data and try to get nodal features as a mesh
include("nodal_features_from_meshing.jl")
grace_factor=10.0
verts′, faces′, freqs, df, (kxs, kys, kzs), Nks = nodal_features_as_mesh(data, bands=(nbandsbelow,nbandsbelow+1); grace_factor=grace_factor)
origin=SVector(minimum(kxs), minimum(kys), minimum(kzs))
widths=SVector(0.5, 0.5, 0.5)
verts′, faces′ = nodal_features_as_mesh_symmetry_extended(df, grace_factor*minimum(df), sgnum, 
    origin=origin, widths=widths, addinv=true)
# verts′ = (verts′ .+ 0.5) .% 1.0 .- 0.5

println(Nks)

# # plotting
# PyPlot.close("all")
# PyPlot.figure()
# plot_trisurf(verts′[:,1], verts′[:,2], verts′[:,3], triangles = faces′ .- 1)
# # for (lab, kv) in highsym_kvs
# #     scatter3D(kv[1], kv[2], kv[3], "o")
# # end
# # xlim(extrema(kxs)), ylim(extrema(kys)), zlim(extrema(kzs))
# xlim((-0.5, 0.5)), ylim((-0.5, 0.5)), zlim((-0.5, 0.5))
# xlabel(L"k_x"), ylabel(L"k_y"), zlabel(L"k_z")
# PyPlot.display_figs()

# using GLMakie
# c = wignerseitz(Gs)
# fig, ax, p = GLMakie.plot(c)
# kp = irrfbz_path(sgnum, Rs)
# GLMakie.plot!(kp)
# # Convert vertices to Cartesian basis
# for (i, vert) in enumerate(eachrow(verts′))
#     verts′[i, :] = Brillouin.cartesianize(vert, Gs)
# end
# GLMakie.mesh!(ax, verts, faces′, color=verts[:, 1])
# display(fig)

verts = verts′
fig = GLMakie.Figure()
ax = Axis3(fig[1, 1], azimuth=-0.2pi, aspect=:data, perspectiveness=0.2)
GLMakie.mesh!(ax, verts, faces′, color=verts[:, 1])
display(fig)

# PyPlot.figure()
# #pcolor(kys, kzs, df[end,:,:], cmap="hot")
# Ncnts = 50
# CLs=contourf(kys, kzs, df[div(Nks[1],2)+1,:,:]', Ncnts, vmin=0.0, vmax=0.055, cmap="inferno")
# colorbar()
# PyPlot.contour(kys, kzs, df[div(Nks[1],2)+1,:,:]', Ncnts, linestyles="solid", colors="white", linewidths=.5)
# #gca().axis("equal")
# #xlim((0.1,0.15)), ylim((0.3, 0.7))
# xlabel(L"k_y"), ylabel(L"k_z")

# PyPlot.figure()
# #pcolor(kys, kzs, df[end,:,:], cmap="hot")
# for bb in 2:3
#     surf(kys, kzs, reshape(freqs[:,bb], Nks)[div(Nks[1],2)+1,:,:]')
# end
# xlabel(L"k_y"), ylabel(L"k_z")
# xlim((0.1,0.15)), ylim((0.3,0.7))
# gca().axis("auto")
# 

# traces, layout = plot_nodal_lines_in_bz(verts′, faces′, sgnum, Rs)
# PlotlyJS.plot(trace, layout=layout)
