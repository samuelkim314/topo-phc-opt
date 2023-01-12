using Crystalline, DelimitedFiles, PlotlyJS, Brillouin, Interpolations
using MPBUtils
include("bz_utils.jl")
include("nodal_features_from_meshing.jl")

sgnum    = 13
subdir = "sg$sgnum-dim$D-n4.00/opt-local-4-579-best"
calcname = "$subdir/3944-bzsymreduced"
subdir = "sg$sgnum-dim$D-n4.00/opt-local-7-1242-best"
calcname = "$subdir/3918-bzsymreduced"
# subdir = "sg$sgnum-dim$D-n4.00/opt-local-9-1916-best"
# calcname = "$subdir/4417-bzsymreduced"

# subdir="sg13-dim3-n4.00/opt-local-4-579-best-3944/local-eps14-best"
# calcname="$subdir/2044-bzsymreduced"
# subdir="sg13-dim3-n4.00/opt-local-4-579-best-3944/local-eps12-best"
# calcname="$subdir/324-bzsymreduced"
# subdir="sg13-dim3-n4.00/opt-local-4-579-best-3944/local-eps10-best"
# calcname="$subdir/3059-bzsymreduced"

outname  = "output/$calcname-dispersion.out"
data = readdlm(outname, ',', Float64)
Rs, flat, isoval, epsin, epsout, _ = lattice_from_mpbparams("input/$calcname-temp.sh") # Rs, flat, isoval, epsin, epsout, kvecs
Gs = reciprocalbasis(Rs)

# Extract data
freqs, (kxs, kys, kzs), (Nkx, Nky, Nkz), (Δx, Δy, Δz) = unpack_mpb_data(data)
freqs = reshape(freqs, (Nkx, Nky, Nkz, size(freqs, 2)))

# Calculate coordinate of nodal line in each ky slice
coords = zeros(Nky, 3)
for i in eachindex(kys)
    freqxz = abs.(freqs[:, i, :, 3] - freqs[:, i, :, 2])    # Band difference
    # Calculate contour around minimum and centroid of contour polygon
    centroidi = calc2dcentroid(freqxz, kxs, kzs, threshold=2.0)
    coords[i, :] = [centroidi[1], kys[i], centroidi[2]]
end

# Interpolate bands to get nodal line frequency in each ky slice
freqnodal = []
for i in eachindex(kys)
    itp2 = LinearInterpolation((kxs, kys), freqs[:, i, :, 2])
    itp3 = LinearInterpolation((kxs, kys), freqs[:, i, :, 3])
    kxi, kzi = coords[i, 1], coords[i, 3]
    freqiavg = (itp2[kxi, kzi] + itp3[kxi, kzi]) / 2
    append!(freqnodal, freqiavg)
end

nodalfreqs = extrema(freqnodal)
println("Min/max nodal line frequency: $nodalfreqs") 
println("Nodal line variation: ", 2*(nodalfreqs[2] - nodalfreqs[1])/(nodalfreqs[2] + nodalfreqs[1]))

p = PlotlyJS.plot(PlotlyJS.scatter(;x=kys, y=freqnodal, mode="line"),
    Layout(;xaxis_title="k_y", yaxis_title="Frequency", xaxis_showgrid=false))
# PlotlyJS.savefig(p, "fullbz/nodalfreq.png")