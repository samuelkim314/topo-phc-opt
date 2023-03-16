using Crystalline, PyPlot, DelimitedFiles, LinearAlgebra, Interpolations
using MPBUtils
isdefined(Main, :KPaths) || include((@__DIR__)*"/../kpaths.jl")
using Main.KPaths

sgnum    = 68
D        = 3
res      = 16
calcidx  = 15
Nk       = 51
calcname = "dim$(D)-sg$(sgnum)-bzfull_$(calcidx)_Nk$(Nk)-res$(res)"
#calcname = "dim$(D)-sg$(sgnum)-bzsymreduced_$(calcidx)_Nk$(Nk)-res$(res)"
outname  = (@__DIR__)*"/../../output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
Rs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams((@__DIR__)*"/../../input/$calcname.sh") # Rs, flat, isoval, epsin, epsout, kvecs
Gs = reciprocalbasis(Rs)
#plot(flat, Rs; isoval=isoval)

# extract data and try to get nodal features as a mesh
include("nodal_features_from_meshing.jl")
verts′, faces′, freqs, df, (kxs, kys, kzs), Nks = nodal_features_as_mesh(data, bands=(2,3))

# high-sym k-points in 68 (primitive) used in symmetry vector; hardcoded
highsym_kvs = Dict("Γ"=>[0.0, 0.0, 0.0], "T"=>[0.5, 0.5, 0.5], "Y"=>[0.5, 0.5, 0.0],
                   "Z"=>[0.0, 0.0, 0.5], "R"=>[0.0, 0.5, 0.5], "S"=>[0.0, 0.5, 0.0])


# plotting
PyPlot.close("all")
PyPlot.figure()
plot_trisurf(verts′[:,1], verts′[:,2], verts′[:,3], triangles = faces′ .- 1)
for (lab, kv) in highsym_kvs
    scatter3D(kv[1], kv[2], kv[3], "o")
end
xlim(extrema(kxs)), ylim(extrema(kys)), zlim(extrema(kzs))
#xlim((0.49, 0.51)), ylim((0.1,0.15)), zlim((0.3,0.7))
xlabel(L"k_x"), ylabel(L"k_y"), zlabel(L"k_z")

PyPlot.figure()
#pcolor(kys, kzs, df[end,:,:], cmap="hot")
Ncnts = 50
CLs=contourf(kys, kzs, df[div(Nks[1],2)+1,:,:]', Ncnts, vmin=0.0, vmax=0.055, cmap="inferno")
colorbar()
contour(kys, kzs, df[div(Nks[1],2)+1,:,:]', Ncnts, linestyles="solid", colors="white", linewidths=.5)
#gca().axis("equal")
#xlim((0.1,0.15)), ylim((0.3, 0.7))
xlabel(L"k_y"), ylabel(L"k_z")

PyPlot.figure()
#pcolor(kys, kzs, df[end,:,:], cmap="hot")
for bb in 2:3
    surf(kys, kzs, reshape(freqs[:,bb], Nks)[div(Nks[1],2)+1,:,:]')
end
xlabel(L"k_y"), ylabel(L"k_z")
xlim((0.1,0.15)), ylim((0.3,0.7))
gca().axis("square")