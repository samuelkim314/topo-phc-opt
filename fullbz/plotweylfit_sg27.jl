using Crystalline, DelimitedFiles, LinearAlgebra, PlotlyJS, Contour, Brillouin, Interpolations
using MPBUtils
isdefined(Main, :KPaths) || include((@__DIR__)*"/../kpaths.jl")
using Main.KPaths

sgnum    = 27
D        = 3
res      = 32
calcidx  = 15
Nk       = 51
# subdir = "sg$sgnum-dim$D-n4.00/ff-direct/rand-local-6-4724-best"
# calcname = "$subdir/3160-bzsymreduced"
# subdir = "sg$sgnum-dim$D-n4.00/ff-direct/rand-local-7-484-best"
# calcname = "$subdir/4969-bzsymreduced"
subdir = "sg$sgnum-dim$D-n4.00/ff-direct/local-8-434-best"
calcname = "$subdir/4823-bzsymreduced"
outname  = "output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
Rs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("output/$calcname-temp.sh") # Rs, flat, isoval, epsin, epsout, kvecs
Gs = reciprocalbasis(Rs)

# Extract data
freqs = sort!(data[:,6:end], dims=2)
kvs = collect(eachrow(data[:,2:4]))
kxs, kys, kzs = sort(unique(getindex.(kvs, 1))), sort(unique(getindex.(kvs, 2))), sort(unique(getindex.(kvs, 3)))
Δx, Δy, Δz = -(-)(extrema(kxs)...), -(-)(extrema(kys)...), -(-)(extrema(kzs)...)
Nkx, Nky, Nkz = length(kxs), length(kys), length(kzs)
freqs = reshape(freqs, (Nkx, Nky, Nkz, size(freqs, 2)))

ikzweyl = Nkz

# Centroid of a small contour to get the location of extrema
freqxy = abs.(freqs[:, :, end, 5] - freqs[:, :, end, 4])
# Calculate contour lines
minline = lines(Contour.contour(kxs, kys, freqxy, 1.5*minimum(freqxy)))
# minline = lines(Contour.contour(kxs, kzs, freqxz, 5.0*minimum(min23)))
contourcoordinates = Contour.coordinates(minline[1])
contourcoordinates2 = [[x, y] for (x, y) in zip(contourcoordinates...)]

cx = 0
cy = 0
a = 0
npoints = size(contourcoordinates2, 1)
for j in 1:npoints-1
    global cx, cy, a
    xi = contourcoordinates2[j][1]
    yi = contourcoordinates2[j][2]
    xi1 = contourcoordinates2[j%npoints+1][1]
    yi1 = contourcoordinates2[j%npoints+1][2]
    aj = xi*yi1 - xi1*yi
    cx += (xi+xi1)*aj
    cy += (yi+yi1)*aj
    a += aj
end
a /= 2
cx /= (6*a)
cy /= (6*a)

centroid = [cx, cy, kzs[end]]

# Apply symmetry operators
coords_list = [centroid]
pg = pointgroup(primitivize(spacegroup(sgnum, Val(3))))
for op in pg
    global coords
    rotation(op) == I && continue
    # coords = vcat(coords, (rotation(op)'\coords')')
    newcoord = rotation(op)'\centroid
    push!(coords_list, newcoord)
end

# Transform to Cartesian basis
coords_list = [Brillouin.cartesianize(coord, Gs) for coord in coords_list]
coords = hcat(coords_list...)'

c = wignerseitz(Gs)
p = Brillouin.plot(c)
BZtrace = p.plot.data
layout = p.plot.layout

nodaltraces = [PlotlyJS.scatter(x=coords[:, 1], y=coords[:, 2], z=coords[:, 3], 
    type="scatter3d", mode="markers", marker_color="green", marker_size=5)]

# append!(nodaltraces, [PlotlyJS.scatter(x=k_interp[:, 1], y=k_interp[:, 2], z=k_interp[:, 3], type="scatter3d")])
# PlotlyJS.plot([PlotlyJS.scatter(x=coords[:, 1], y=coords[:, 2], z=coords[:, 3], type="scatter3d", mode="lines", marker_color="black")],
#     Layout(scene_xaxis_range=[-0.5, 0.5], scene_yaxis_range=[-0.5, 0.5], scene_zaxis_range=[-0.5, 0.5]))
# PlotlyJS.plot(nodaltraces,
#     Layout(scene_xaxis_range=[-0.5, 0.5], scene_yaxis_range=[-0.5, 0.5], scene_zaxis_range=[-0.5, 0.5]))
p = PlotlyJS.plot(append!(nodaltraces, BZtrace),
    p.plot.layout)
# PlotlyJS.savefig(p, "fullbz/sg27weyl.png")