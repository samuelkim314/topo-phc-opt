using Crystalline, DelimitedFiles, PlotlyJS
using MPBUtils
include("nodal_features_from_meshing.jl")

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

# extract data and try to get nodal features as a mesh
grace_factor=5.0
verts′, faces′, freqs, df, (kxs, kys, kzs), Nks = nodal_features_as_mesh(data, bands=(2,3); grace_factor=grace_factor)


verts′ = sortslices(verts′, dims=1, by=x->x[2]) # Sort rows based on kys value
freqs = sort!(data[:,6:end], dims=2)
freqs = reshape(freqs, (Nks[1], Nks[2], Nks[3], size(freqs, 2)))
nodalfreqs = []
nodalfreqs2 = []
vertkys = verts′[:, 2]
for vert in eachrow(verts′)
    ix = argmin(abs.(kxs .- vert[1]))
    iy = argmin(abs.(kys .- vert[2]))
    iz = argmin(abs.(kzs .- vert[3]))
    append!(nodalfreqs, freqs[ix, iy, iz, 2])
    append!(nodalfreqs2, freqs[ix, iy, iz, 3])
end
# t1 = PlotlyJS.scatter(;x=vertkys, y=nodalfreqs, mode="line", name="Band 2")
# t2 = PlotlyJS.scatter(;x=vertkys, y=nodalfreqs2, mode="line", name="Band 3")
# p = PlotlyJS.plot([t1, t2])
# PlotlyJS.savefig(p, "fullbz/nodalbandsfreq.png")

p = PlotlyJS.plot(PlotlyJS.scatter(;x=vertkys, y=(nodalfreqs2 .+ nodalfreqs)./2, mode="line"),
    Layout(;xaxis_title="k_y", yaxis_title="Frequency", xaxis_showgrid=false))
PlotlyJS.savefig(p, "fullbz/nodalbandsfreqavg.png")