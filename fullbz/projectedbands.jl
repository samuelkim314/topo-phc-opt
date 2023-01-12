using PyPlot


subdir = "sg$sgnum-dim$D-n4.00/opt-local-9-1916-best"
calcname = "$subdir/4417-bzsymreduced"
subdir = "sg75-dim3-n4.00/chern1-ff-direct/local-8-4919-best"
calcname = "$subdir/4958-bzsymreduced"
subdir = "sg168-dim3-n4.00/chern1-ff-direct/local-5-4611-best"
calcname = "$subdir/2355-bzsymreduced"
outname  = "output/$calcname-dispersion.out"

bands = Vector{Float64}[]

data = readdlm(outname, ',', Float64)
freqs = sort!(data[:,6:end], dims=2)
nbands = size(freqs, 2)
kvs = collect(eachrow(data[:,2:4]))
kxs, kys, kzs = sort(unique(getindex.(kvs, 1))), sort(unique(getindex.(kvs, 2))), sort(unique(getindex.(kvs, 3)))
Δx, Δy, Δz = -(-)(extrema(kxs)...), -(-)(extrema(kys)...), -(-)(extrema(kzs)...)
Nkx, Nky, Nkz = length(kxs), length(kys), length(kzs)
freqs = reshape(freqs, (Nkx, Nky, Nkz, nbands))

println((Nkx, Nky, Nkz, nbands))

for ix in 1:Nkx
    for iy in 1:Nky
        for ib in 1:nbands
            push!(bands, freqs[ix, iy, :, ib])
        end
    end
end
println(size(bands))

PyPlot.figure()
for (i, line) in enumerate(bands)
    icolor = i % nbands
    PyPlot.plot(line, color="C$icolor")
end
PyPlot.display_figs()


# PyPlot.savefig("fullbz/hinge-projected.png")

