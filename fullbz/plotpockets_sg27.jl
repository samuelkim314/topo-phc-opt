using Crystalline, PyPlot, DelimitedFiles, Brillouin
using MPBUtils
include("nodal_features_from_meshing.jl")
include("bz_utils.jl")

makie = true    # Otherwise reverts to PyPlot, but this doesn't plot the BZ
if makie
    using GLMakie
    # using CairoMakie
    # GLMakie.activate!(ssao=true)
end


# Input parameters
sgnum    = 27
D        = 3
# subdir = "sg$sgnum-dim$D-n4.00/ff-direct/rand-local-6-4724-best"
# id = 3160
# subdir = "sg$sgnum-dim$D-n4.00/ff-direct/rand-local-7-484-best"
# id = 4969
subdir = "sg$sgnum-dim$D-n4.00/ff-direct/local-8-434-best"
id = 4823
calcname = "$subdir/$id-bzsymreduced"
nbandsbelow = 4
outname  = "output/$calcname-dispersion.out"

# Read the data
data = readdlm(outname, ',', Float64)
Rs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("output/$subdir/$id.sh") # Rs, flat, isoval, epsin, epsout, kvecs
Gs = reciprocalbasis(Rs)

# Centroid of a small contour to get the location of extrema
freqs, (kxs, kys, kzs), (Nkx, Nky, Nkz), (Δx, Δy, Δz) = unpack_mpb_data(data)
freqs = reshape(freqs, Nkx, Nky, Nkz, size(freqs, 2))

ikzweyl = Nkz   # Weyl points are at kz=pi plane (assume BZ simulation extends to kz=pi)
freqxy = abs.(freqs[:, :, ikzweyl, nbandsbelow+1] - freqs[:, :, ikzweyl, nbandsbelow])
# Calculate centroid from the contour
centroid = calc2dcentroid(freqxy, kxs, kys)     # centroid in xy plane
append!(centroid, kzs[ikzweyl])     # Convert to 3D coordinate
println("Centroid: $centroid")

# Apply symmetry operators and inversion operator
coords = applysymmetry(centroid, sgnum, true)
# Convert Weyl point coordinates to Cartesian basis
for (i, coord) in enumerate(eachrow(coords))
    coords[i, :] = Brillouin.cartesianize(coord, Gs)
end

# Interpolate bands at k_3=0.5 plane to get Weyl point frequency
weylfreq = getweylfreq(centroid, freqs, kxs, kys, kzs, nbandsbelow)

# Calculate mesh for Fermi pockets
verts, faces = isosurfacemesh(data, weylfreq, nbandsbelow, true, true, invertfaces=true)

c = wignerseitz(Gs)
kp = irrfbz_path(sgnum, Rs)

if makie
    fig, ax, p = GLMakie.plot(c)
    # GLMakie.mesh!(ax, verts, faces, color=verts[:, 1])
    GLMakie.mesh!(ax, verts, faces, shading=true, rasterize=5, color=:yellowgreen, transparency=false,
        shininess=32f0)
    GLMakie.plot!(kp, textkws=(; strokewidth=0))

    # plot_trisurf(verts[:,1], verts[:,2], verts[:,3], triangles = faces .- 1)

    GLMakie.scatter!(ax, coords[:, 1], coords[:, 2], coords[:, 3]; markersize=15)
    # ax.viewmode = :fit
    ax.azimuth = 0.15pi
    display(fig)
    # CairoMakie.save("figures/sg27-4823-fermi.pdf", fig)
else
    PyPlot.close("all")
    PyPlot.figure()
    plot_trisurf(verts[:,1], verts[:,2], verts[:,3], triangles = faces .- 1)
    PyPlot.scatter3D(coords[:, 1], coords[:, 2], coords[:, 3], s=10)
    # xlim(extrema(kxs)), ylim(extrema(kys)), zlim(extrema(kzs))
    # xlim((-0.5, 0.5)), ylim((-0.5, 0.5)), zlim((-0.5, 0.5))
    xlim((-pi, pi)), ylim((-pi, pi)), zlim((-pi, pi))
    xlabel(L"k_1"), ylabel(L"k_2"), zlabel(L"k_3")

    PyPlot.display_figs()

end

