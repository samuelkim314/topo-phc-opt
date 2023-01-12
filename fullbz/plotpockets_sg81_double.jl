#=
Plot Fermi pockets for SG 81, for the case where there are 8 Weyl points
(4 at kz=0 and 4 at kz=pi)
=#
using Crystalline, PyPlot, DelimitedFiles, GLMakie, Brillouin
using MPBUtils
include("nodal_features_from_meshing.jl")
include("bz_utils.jl")

sgnum    = 81
D        = 3
subcalcname = "sg$sgnum-dim$D-n4.00/gaps-ff-direct/rand-local-10-4822-best/4358"
nbandsbelow = 7
subcalcname = "sg$sgnum-dim$D-n4.00/gaps-ff-direct/rand-local-9-3098-best/4990"
nbandsbelow = 7

# Load data
calcname = "$subcalcname-bzsymreduced"
outname  = "output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
Rs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("output/$subcalcname.sh") # Rs, flat, isoval, epsin, epsout, kvecs
Gs = reciprocalbasis(Rs)

# Centroid of a small contour to get the location of extrema
freqs, (kxs, kys, kzs), (Nkx, Nky, Nkz), (Δx, Δy, Δz) = unpack_mpb_data(data)
freqs = reshape(freqs, Nkx, Nky, Nkz, size(freqs, 2))

ikzweyl = 1     # Weyl point at kz=0 plane, assume BZ simulation started at kz=0)
freqxy = abs.(freqs[:, :, ikzweyl, nbandsbelow+1] - freqs[:, :, ikzweyl, nbandsbelow])
# Calculate contour and centroid to get Weyl point locaiton
centroid = calc2dcentroid(freqxy, kxs, kys)     # centroid in xy plane
append!(centroid, kzs[ikzweyl])     # Convert to 3D coordinate
println("Centroid 1: $centroid")

# Interpolate bands to get Weyl frequency
weylfreq1 = getweylfreq(centroid, freqs, kxs, kys, kzs, nbandsbelow)

# Apply symmetry operators
coords1 = applysymmetry(centroid, sgnum, true)
# Convert Weyl point coordinates to Cartesian basis
for (i, coord) in enumerate(eachrow(coords1))
    coords1[i, :] = Brillouin.cartesianize(coord, Gs)
end


ikzweyl = Nkz
freqxy = abs.(freqs[:, :, ikzweyl, nbandsbelow+1] - freqs[:, :, ikzweyl, nbandsbelow])
# Calculate centroid from the contour
centroid = calc2dcentroid(freqxy, kxs, kys)     # centroid in xy plane
append!(centroid, kzs[ikzweyl])     # Convert to 3D coordinate
println("Centroid 2: $centroid")

# Calculate second Weyl frequency
weylfreq2 = getweylfreq(centroid, freqs, kxs, kys, kzs, nbandsbelow)

# Apply symmetry operators
coords2 = applysymmetry(centroid, sgnum, true)
# Convert Weyl point coordinates to Cartesian basis
for (i, coord) in enumerate(eachrow(coords2))
    coords2[i, :] = Brillouin.cartesianize(coord, Gs)
end


# Calculate mesh for Fermi pockets
weylfreq = (weylfreq1 + weylfreq2) / 2
verts, faces = isosurfacemesh(data, weylfreq, nbandsbelow, true, true, invertfaces=true)

# Plot
c = wignerseitz(Gs)
fig, ax, p = GLMakie.plot(c)
GLMakie.mesh!(ax, verts, faces, color=verts[:, 1])
kp = irrfbz_path(sgnum, Rs)
GLMakie.plot!(kp)
# display(fig)

# plot_trisurf(verts[:,1], verts[:,2], verts[:,3], triangles = faces .- 1)

GLMakie.scatter!(ax, coords1[:, 1], coords1[:, 2], coords1[:, 3]; markersize=15)
GLMakie.scatter!(ax, coords2[:, 1], coords2[:, 2], coords2[:, 3]; markersize=15)
# ax.aspect = (1, 1, 1)
# ax.aspect = :data
# ax.viewmode = :fit
display(fig)

# PyPlot.scatter3D(coords[:, 1], coords[:, 2], coords[:, 3], s=10)

