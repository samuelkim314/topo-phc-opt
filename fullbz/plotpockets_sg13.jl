using Crystalline, PyPlot, DelimitedFiles, LinearAlgebra, Interpolations, PlotlyJS, Meshing
using MPBUtils
using Main.KPaths
using Crystalline: _mesh_to_cartesian, ReciprocalBasis
include("nodal_features_from_meshing.jl")

sgnum    = 13
D        = 3
subdir = "sg$sgnum-dim$D-n4.00/opt-local-4-579-best"
calcname = "$subdir/3944-bzsymreduced"
outname  = "output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
Rs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("input/$calcname-temp.sh") # Rs, flat, isoval, epsin, epsout, kvecs
Gs = reciprocalbasis(Rs)

# plotting
PyPlot.close("all")
PyPlot.figure()

# isoval = 0.3287304
isoval = 0.317
# isoval = 0.3154623

for isoval in range(0.3154623, 0.3287304, length=3)
    freqs, (kxs, kys, kzs), (Nkx, Nky, Nkz), (Δx, Δy, Δz) = unpack_mpb_data(data)
    origin=SVector(minimum(kxs), minimum(kys), minimum(kzs))
    widths=SVector(Δx, Δy, Δz)
    # verts, faces = Meshing.isosurface(freqs4, MarchingTetrahedra(iso=isoval, eps=1e-5), origin=origin, widths=widths)
    # Gstemp=ReciprocalBasis{3}([1.,0,0], [0,1.,0], [0,0,1.])
    # verts, faces = _mesh_to_cartesian(verts, faces, Gstemp)
    freqs2 = reshape(freqs[:,2], Nkx, Nky, Nkz)
    freqs3 = reshape(freqs[:,3], Nkx, Nky, Nkz)
    verts2, faces2 = nodal_features_as_mesh_symmetry_extended(freqs2, isoval, sgnum; origin=origin, widths=widths)
    verts3, faces3 = nodal_features_as_mesh_symmetry_extended(freqs3, isoval, sgnum; origin=origin, widths=widths)
    verts = vcat(verts2, verts3)
    faces = vcat(faces2, faces3.+size(verts2, 1))
    println(size(verts))

    plot_trisurf(verts[:,1], verts[:,2], verts[:,3], triangles = faces .- 1)
end


# xlim(extrema(kxs)), ylim(extrema(kys)), zlim(extrema(kzs))
xlim((-0.5, 0.5)), ylim((-0.5, 0.5)), zlim((-0.5, 0.5))
xlabel(L"k_1"), ylabel(L"k_2"), zlabel(L"k_3")

# PyPlot.display_figs()
PyPlot.savefig("fullbz/figs/sg13pockets.png")
