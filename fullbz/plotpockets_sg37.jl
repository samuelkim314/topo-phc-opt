using Crystalline, PyPlot, DelimitedFiles
using Brillouin
using MPBUtils
include("bz_utils.jl")
include("nodal_features_from_meshing.jl")

makie = true    # Otherwise reverts to PyPlot, but this doesn't plot the BZ
if makie
    using GLMakie
    # using CairoMakie
    # GLMakie.activate!(ssao=true)
end

sgnum    = 37
D        = 3
subcalcname = "sg$sgnum-dim$D-n4.00/ff-direct/rand-local-6-904-best/3730"
# subcalcname = "sg37-dim$D-n4.00/ff-direct/rand-local-9-4279-best/4995"
# subcalcname = "sg$sgnum-dim$D-n4.00/ff-direct/local-2-517-best/2103"
nbandsbelow = 4

calcname = "$subcalcname-bzsymreduced"
outname  = "output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
pRs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("output/$subcalcname.sh") # Rs, flat, isoval, epsin, epsout, kvecs
cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)
Rs = conventionalize(pRs, cntr)
Gs = reciprocalbasis(pRs)

freqs, (kxs, kys, kzs), (Nkx, Nky, Nkz), (Δx, Δy, Δz) = unpack_mpb_data(data)
freqs = reshape(freqs, Nkx, Nky, Nkz, size(freqs, 2))

ikzweyl = Nkz   # Weyl points are at kz=pi plane (assume BZ simulation extends to kz=pi)
freqxy = abs.(freqs[:, :, ikzweyl, nbandsbelow+1] - freqs[:, :, ikzweyl, nbandsbelow])
# Calculate centroid from the contour
centroid = calc2dcentroid(freqxy, kxs, kys)
append!(centroid, kzs[ikzweyl])
println("Centroid: $centroid")
# Interpolate bands at k_3=0.5 plane to get Weyl point frequency
weylfreq = getweylfreq(centroid, freqs, kxs, kys, kzs, nbandsbelow)

# Apply symmetry operators and inversion operator
centroid = Array(Brillouin.reduce_to_wignerseitz(centroid, Gs))
coords = applysymmetry(centroid, sgnum, true)

# Convert Weyl point coordinates to Cartesian basis
for (i, coord) in enumerate(eachrow(coords))
    coords[i, :] = Brillouin.cartesianize(coord, Gs)
end

# Plot Fermi pockets
verts, faces = isosurfacemesh(data, weylfreq, nbandsbelow, true, invertfaces=false)

c = wignerseitz(Gs)
kp = irrfbz_path(sgnum, Rs)

if makie

    fig, ax, p = GLMakie.plot(c)
    # GLMakie.mesh!(ax, verts, faces, color=verts[:, 1])
    GLMakie.mesh!(ax, verts, faces, shading=true, rasterize=5, color=:yellowgreen, transparency=false,
        shininess=32f0)
    # Plot high-symmetry k-path

    GLMakie.plot!(kp, textkws=(; strokewidth=0))
    ax.azimuth = 0.1pi

    # Plot Weyl point coordinates
    GLMakie.scatter!(ax, coords[:, 1], coords[:, 2], coords[:, 3]; markersize=15)
    # ax.azimuth = 0.2pi
    display(fig)
    # CairoMakie.save("figures/sg37-3730-fermi.pdf", fig)
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