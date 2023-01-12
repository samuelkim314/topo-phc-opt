using Crystalline, DelimitedFiles, LinearAlgebra, PlotlyJS, Contour, Brillouin
using MPBUtils
include("bz_utils.jl")

sgnum    = 13
subdir = "sg$sgnum-dim$D-n4.00/opt-local-4-579-best"
calcname = "$subdir/3944-bzsymreduced"
outname  = "output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
Rs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("input/$calcname-temp.sh") # Rs, flat, isoval, epsin, epsout, kvecs
Gs = reciprocalbasis(Rs)
#plot(flat, Rs; isoval=isoval)

# println(irrfbz_path(sgnum, Rs))
# println(Gs)

freqs = sort!(data[:,6:end], dims=2)
kvs = collect(eachrow(data[:,2:4]))
kxs, kys, kzs = sort(unique(getindex.(kvs, 1))), sort(unique(getindex.(kvs, 2))), sort(unique(getindex.(kvs, 3)))
Δx, Δy, Δz = -(-)(extrema(kxs)...), -(-)(extrema(kys)...), -(-)(extrema(kzs)...)
Nkx, Nky, Nkz = length(kxs), length(kys), length(kzs)
freqs = reshape(freqs, (Nkx, Nky, Nkz, size(freqs, 2)))

i = 30
freqxz = abs.(freqs[:, i, :, 3] - freqs[:, i, :, 2])
centroid = calc2dcentroid(freqxz, kxs, kzs, threshold=2.0)

cx = centroid[1]
cz = centroid[2]
println("Centroid: $centroid")

PlotlyJS.plot([PlotlyJS.contour(z=transpose(freqxz), x=kxs, y=kzs),
    PlotlyJS.scatter(x=contourcoordinates[1], y=contourcoordinates[2], mode="line"),
    PlotlyJS.scatter(x=[centroid[1]], y=[centroid[2]], mode="markers"),
    PlotlyJS.scatter(x=[cx], y=[cz], mode="markers")])