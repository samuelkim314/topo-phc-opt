using Crystalline, PyPlot, DelimitedFiles, LinearAlgebra, Brillouin
using MPBUtils
using Crystalline: _mesh_to_cartesian, ReciprocalBasis
using GLMakie
# using CairoMakie
include("nodal_features_from_meshing.jl")
include("bz_utils.jl")

sgnum    = 81
D        = 3
subcalcname = "sg$sgnum-dim$D-n4.00/gaps-ff-direct/local-9-3098-best/4384"
nbandsbelow = 4
# subcalcname = "sg$sgnum-dim$D-n4.00/gaps-ff-direct/local-7-349-best/2840"
# nbandsbelow = 4

calcname = "$subcalcname-bzsymreduced"
outname  = "output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
Rs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("output/$subcalcname.sh") # Rs, flat, isoval, epsin, epsout, kvecs
Gs = reciprocalbasis(Rs)

# Centroid of a small contour to get the location of extrema
freqs, (kxs, kys, kzs), (Nkx, Nky, Nkz), (Δx, Δy, Δz) = unpack_mpb_data(data)
freqs = reshape(freqs, Nkx, Nky, Nkz, size(freqs, 2))

ikzweyl = Nkz
freqxy = abs.(freqs[:, :, ikzweyl, nbandsbelow+1] - freqs[:, :, ikzweyl, nbandsbelow])
# Calculate centroid from the contour
centroid = calc2dcentroid(freqxy, kxs, kys)     # centroid in xy plane
append!(centroid, kzs[ikzweyl])     # Convert to 3D coordinate
println("Centroid: $centroid")

# # Apply symmetry operators, including inversion
coords = applysymmetry(centroid, sgnum, true)

# Convert Weyl point coordinates to Cartesian basis
for (i, coord) in enumerate(eachrow(coords))
    coords[i, :] = Brillouin.cartesianize(coord, Gs)
end

# Interpolate bands at k_3=0.5 plane to get Weyl point frequency
weylfreq = getweylfreq(centroid, freqs, kxs, kys, kzs, nbandsbelow)

# Calculate mesh for Fermi pockets
verts, faces = isosurfacemesh(data, weylfreq, nbandsbelow, true, true)

# fig = GLMakie.Figure()
# ax = Axis3(fig[1, 1], azimuth=-0.2pi, aspect=:data, perspectiveness=0.2)
c = wignerseitz(Gs)
fig, ax, p = GLMakie.plot(c)
if length(verts) > 0
    GLMakie.mesh!(ax, verts, faces, color=verts[:, 1])
end
kp = irrfbz_path(sgnum, Rs)
GLMakie.plot!(kp)
# display(fig)

# plot_trisurf(verts[:,1], verts[:,2], verts[:,3], triangles = faces .- 1)

GLMakie.scatter!(ax, coords[:, 1], coords[:, 2], coords[:, 3]; markersize=15)
# ax.aspect = (1, 1, 1)
# ax.aspect = :data
# ax.viewmode = :fit
ax.azimuth = 0.2pi
display(fig)

# CairoMakie.save("plot.pdf", fig)


# PyPlot.scatter3D(coords[:, 1], coords[:, 2], coords[:, 3], s=10)

