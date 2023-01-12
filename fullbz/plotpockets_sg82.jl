using Crystalline, PyPlot, DelimitedFiles, LinearAlgebra, Meshing, GLMakie, Brillouin
using MPBUtils
# using CairoMakie
using Crystalline: _mesh_to_cartesian, ReciprocalBasis
include("nodal_features_from_meshing.jl")
include("bz_utils.jl")

sgnum    = 82
D        = 3
# subdir = "sg82-dim$D-n4.00/ff-direct/local-1-4286-best"
# calcname = "$subdir/4997-bzsymreduced"
# nbandsbelow = 6
# subdir = "sg$sgnum-dim$D-n4.00/ff-direct/rand-local-4-3886-best"
# calcname = "$subdir/4803-bzsymreduced"
# nbandsbelow = 6
# subdir = "sg$sgnum-dim$D-n4.00/ff-direct/local-2-4911-best"
# calcname = "$subdir/3136-bzsymreduced"
# nbandsbelow = 6
# subdir = "sg$sgnum-dim$D-n4.00/ff-direct/local-4-3886-best"
# calcname = "$subdir/1690-bzsymreduced"
# nbandsbelow = 6
subdir = "sg$sgnum-dim$D-n4.00/ff-direct/rand-local-8-4882-best"
calcname = "$subdir/4920-bzsymreduced"
nbandsbelow = 7
outname  = "output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
pRs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("output/$subdir/4920.sh") # Rs, flat, isoval, epsin, epsout, kvecs
cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)
Rs = conventionalize(pRs, cntr)
Gs = reciprocalbasis(pRs)

freqs, (kxs, kys, kzs), (Nkx, Nky, Nkz), (Δx, Δy, Δz) = unpack_mpb_data(data)
println(size(freqs))

# Centroid of mesh polyhedron
grace_factor=5.0
verts, faces, freqs, df, (kxs, kys, kzs), Nks = nodal_features_as_mesh(data, bands=(nbandsbelow,nbandsbelow+1); grace_factor=grace_factor)

# Sort through verts/faces to pick out one Weyl point (or else centroid calculation is nonsense)
faces = faces[verts[faces[:, 3], 3] .> 0.0, :]

# Calculate centroid of polyhedron
centroid = calc3dcentroid(verts, faces)

println(centroid)
println(Brillouin.cartesianize(centroid, Gs))

# Get Weyl point frequency
freqs = reshape(freqs, Nkx, Nky, Nkz, size(freqs, 2))
weylfreq = getweylfreq(centroid, freqs, kxs, kys, kzs, nbandsbelow)

# Apply symmetry operators
coords = applysymmetry(centroid, sgnum, true)

# Convert Weyl point coordinates to Cartesian basis
for (i, coord) in enumerate(eachrow(coords))
    coords[i, :] = Brillouin.cartesianize(coord, Gs)
end

# extract data and try to get nodal features as a mesh
include("nodal_features_from_meshing.jl")
# Plot Fermi pockets
isoval = weylfreq + 0.001
origin=SVector(minimum(kxs), minimum(kys), minimum(kzs))
widths=SVector(Δx, Δy, Δz)
# verts, faces = Meshing.isosurface(freqs4, MarchingTetrahedra(iso=isoval, eps=1e-5), origin=origin, widths=widths)
# Gstemp=ReciprocalBasis{3}([1.,0,0], [0,1.,0], [0,0,1.])
# verts, faces = _mesh_to_cartesian(verts, faces, Gstemp)
freqs, (kxs, kys, kzs), (Nkx, Nky, Nkz), (Δx, Δy, Δz) = unpack_mpb_data(data)
freqs4 = reshape(freqs[:,nbandsbelow], Nkx, Nky, Nkz)
freqs5 = reshape(freqs[:,nbandsbelow+1], Nkx, Nky, Nkz)
verts4, faces4 = nodal_features_as_mesh_symmetry_extended(freqs4, isoval, sgnum; origin=origin, widths=widths, addinv=true)
verts5, faces5 = nodal_features_as_mesh_symmetry_extended(freqs5, isoval, sgnum; origin=origin, widths=widths, addinv=true)
verts = vcat(verts4, verts5)
faces = vcat(faces4, faces5.+size(verts4, 1))

# Convert vertices to Cartesian basis
for (i, vert) in enumerate(eachrow(verts))
    verts[i, :] = Brillouin.cartesianize(vert, Gs)
end

# fig = GLMakie.Figure()
# ax = Axis3(fig[1, 1], azimuth=-0.2pi, aspect=:data, perspectiveness=0.2)
c = wignerseitz(Gs)
fig, ax, p = GLMakie.plot(c)
if length(verts) > 0
    GLMakie.mesh!(ax, verts, faces, color=verts[:, 1])
end
kp = irrfbz_path(sgnum, Rs)
GLMakie.plot!(kp)
# Plot Weyl points
# println(size(coords))
GLMakie.scatter!(ax, coords[:, 1], coords[:, 2], coords[:, 3]; markersize=15)
# ax.azimuth = -0.1pi
display(fig)
# CairoMakie.save("output/$subdir/bz-3136.pdf", fig)


# plot_trisurf(verts[:,1], verts[:,2], verts[:,3], triangles = faces .- 1)
# PyPlot.scatter3D(coords[:, 1], coords[:, 2], coords[:, 3], s=10)

