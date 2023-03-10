using Crystalline, DelimitedFiles, LinearAlgebra, PlotlyJS, Contour, Brillouin, Interpolations
using MPBUtils
isdefined(Main, :KPaths) || include((@__DIR__)*"/../kpaths.jl")
using Main.KPaths

sgnum    = 27
D        = 3
res      = 32
subdir = "sg27-dim$D-n4.00/ff-direct/local-6-4724-best"
calcname = "$subdir/3160-bzsymreduced"
subdir = "sg27-dim$D-n4.00/ff-direct/local-7-484-best"
calcname = "$subdir/4969-bzsymreduced"
# subdir = "sg37-dim$D-n4.00/ff-direct/local-6-904-best"
# calcname = "$subdir/3730-bzsymreduced"
# subdir = "sg37-dim$D-n4.00/ff-direct/local-9-4279-best"
# calcname = "$subdir/4995-bzsymreduced"
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
println((cx, cy))


itp4 = LinearInterpolation((kxs, kys), freqs[:, :, end, 4])
itp5 = LinearInterpolation((kxs, kys), freqs[:, :, end, 5])
weylfreq4 = itp4(cx, cy)
weylfreq5 = itp5(cx, cy)
freqiavg = (weylfreq4 + weylfreq5) / 2

println((weylfreq4, weylfreq5, freqiavg))
# p = PlotlyJS.plot([heatmap(x=kxs, y=kys, z=freqxy'),
#     PlotlyJS.scatter(x=contourcoordinates[1], y=contourcoordinates[2], mode="line"),
#     PlotlyJS.scatter(x=cx, y=cy, mode="markers", marker_size=40,)])
# PlotlyJS.savefig(p, "fullbz/sg37weylslice.png")

# p = PlotlyJS.plot([heatmap(x=kxs, y=kys, z=freqs[:, :, end, 4]'),
#     PlotlyJS.scatter(x=contourcoordinates[1], y=contourcoordinates[2], mode="line"),
#     PlotlyJS.scatter(x=cx, y=cy, mode="markers", marker_size=40,)])
# PlotlyJS.savefig(p, "fullbz/sg37band4slice.png")

# p = PlotlyJS.plot([heatmap(x=kxs, y=kys, z=freqs[:, :, end, 5]'),
#     PlotlyJS.scatter(x=contourcoordinates[1], y=contourcoordinates[2], mode="line"),
#     PlotlyJS.scatter(x=cx, y=cy, mode="markers", marker_size=40,)])
# PlotlyJS.savefig(p, "fullbz/sg37band5slice.png")
