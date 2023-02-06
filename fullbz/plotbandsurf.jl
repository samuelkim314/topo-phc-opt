include("../src/bands.jl")
using PlotlyJS


Nk = 50
kxs = range(0, 0.5, length=Nk)
kys = range(0.0, 1/3,  length=Nk)
kzs = range(0.0, 0.5,  length=Nk)
kvs = [[kx, ky, kz] for kx in kxs, ky in kys, kz in kzs]

println(size(kvs))

topofullpath = "thomas-candidates/dim3-sg147/1to6-local-14674-best/2962-bzsymreduced"

lines = Vector{Float64}[]
freqs = getbands(topofullpath)
freqs = reshape(freqs, (50, 50, 50, 10))

band6 = freqs[:, :, end, 6]
band7 = freqs[:, :, end, 7]

# PyPlot.figure()
# PyPlot.surf(kxs, kys, band6)
# PyPlot.surf(kxs, kys, band7)
# PyPlot.savefig("fullbz/hinge-surface.png")
# PyPlot.display_figs()

PlotlyJS.plot([surface(z=band6, x=kxs, y=kys),
    surface(z=band7, x=kxs, y=kys)])

