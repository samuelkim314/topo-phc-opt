using Crystalline, PyPlot, DelimitedFiles, LinearAlgebra, Interpolations, PlotlyJS, Brillouin
using MPBUtils
using Main.KPaths

sgnum    = 13
D        = 3
res      = 16
calcidx  = 15
Nk       = 51
subdir = "sg$sgnum-dim$D-n4.00/opt-local-4-579-best"
calcname = "$subdir/3944-bzsymreduced"
subdir = "sg$sgnum-dim$D-n4.00/opt-local-7-1242-best"
calcname = "$subdir/3918-bzsymreduced"
# subdir = "sg$sgnum-dim$D-n4.00/opt-local-9-1916-best"
# calcname = "$subdir/4417-bzsymreduced"
subdir="sg13-dim3-n4.00/opt-local-4-579-best-3944/local-eps14-best"
calcname="$subdir/2044-bzsymreduced"
# subdir="sg13-dim3-n4.00/opt-local-4-579-best-3944/local-eps12-best"
# calcname="$subdir/324-bzsymreduced"
# subdir="sg13-dim3-n4.00/opt-local-4-579-best-3944/local-eps10-best"
# calcname="$subdir/3059-bzsymreduced"

outname  = "output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
Rs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("output/$calcname-temp.sh") # Rs, flat, isoval, epsin, epsout, kvecs
Gs = reciprocalbasis(Rs)
#plot(flat, Rs; isoval=isoval)

# println(irrfbz_path(sgnum, Rs))
# println(Gs)

# extract data and try to get nodal features as a mesh
include("nodal_features_from_meshing.jl")
grace_factor=10.0
verts′, faces′, freqs, df, (kxs, kys, kzs), Nks = nodal_features_as_mesh(data, bands=(2,3); grace_factor=grace_factor)
origin=SVector(minimum(kxs), minimum(kys), minimum(kzs))
widths=SVector(1.0, 0.5, 0.5)
verts′, faces′ = nodal_features_as_mesh_symmetry_extended(df, grace_factor*minimum(df), sgnum, 
    origin=origin, widths=widths)

# high-sym k-points in 68 (primitive) used in symmetry vector; hardcoded
highsym_kvs = Dict("Γ"=>[0.0, 0.0, 0.0], "D"=>[0.0, 0.5, 0.5], "B"=>[0.0, 0.0, 0.5],
                   "Z"=>[0.0, 0.5, 0.0], "A"=>[-0.5, 0.0, 0.5], "E"=>[-0.5, 0.5, 0.5],
                   "C₂"=>[-0.5, 0.5, 0.0], "Y₂"=>[-0.5, 0.0, 0.0])


# plotting
PyPlot.close("all")
PyPlot.figure()
plot_trisurf(verts′[:,1], verts′[:,2], verts′[:,3], triangles = faces′ .- 1)
# for (lab, kv) in highsym_kvs
#     scatter3D(kv[1], kv[2], kv[3], "o")
# end
# xlim(extrema(kxs)), ylim(extrema(kys)), zlim(extrema(kzs))
xlim((-0.5, 0.5)), ylim((-0.5, 0.5)), zlim((-0.5, 0.5))
#xlim((0.49, 0.51)), ylim((0.1,0.15)), zlim((0.3,0.7))
xlabel(L"k_x"), ylabel(L"k_y"), zlabel(L"k_z")

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
PyPlot.display_figs()

# traces, layout = plot_nodal_lines_in_bz(verts′, faces′, sgnum, Rs)
# PlotlyJS.plot(trace, layout=layout)
