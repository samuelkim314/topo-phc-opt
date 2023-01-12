using DelimitedFiles
using MPBUtils, Crystalline, Brillouin

sgnum    = 82
D        = 3
# calcname = "sg13-dim$D-n4.00/opt-local-9-1916-best/4417"
# calcname = "sg$sgnum-dim$D-n4.00/ff-direct/local-6-4724-best/3160"
# calcname = "sg82-dim$D-n4.00/ff-direct/local-1-4286-best/4997"
# calcname = "sg82-dim$D-n4.00/ff-direct/local-2-4911-best/3136"
# calcname = "sg82-dim$D-n4.00/ff-direct/rand-local-1-4286-best/4997"
# calcname = "sg82-dim$D-n4.00/ff-direct/rand-local-4-3886-best/4803"
# calcname = "sg82-dim$D-n4.00/ff-direct/rand-local-8-4882-best/4920"
# calcname = "sg81-dim$D-n4.00/gaps-ff-direct/local-9-3098-best/4990"
# calcname = "sg81-dim$D-n4.00/gaps-ff-direct/local-10-4822-best/4358"
# calcname = "sg81-dim$D-n4.00/gaps-ff-direct/local-9-3098-best/4384"
# outname  = "output/$calcname-dispersion.out"
# data = readdlm(outname, ',', Float64)
# data = data[:, 1:end .!= 5]
# data = data[:, 2:end]
# writedlm("fullbz/DOS-calculation/sg$sgnum-8-4882-4920-band.txt", data, ' ')

outname  = "output/$calcname-bzsymreduced-dispersion.out"
data = readdlm(outname, ',', Float64)
data = data[:, 1:end .!= 5]     # Remove column corresponding to k-point magnitude
data = data[:, 2:end]           # Remove column for k-point index
writedlm("fullbz/DOS-calculation/data/sg$sgnum-2-4911-3136-freq_Tr.txt", data, ' ')

pRs, _, _, _, _, _ = lattice_from_mpbparams("output/$calcname.sh")
cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)
Rs = conventionalize(pRs, cntr)
Gs = reciprocalbasis(pRs)
println(Gs)
writedlm("fullbz/DOS-calculation/data/sg$sgnum-2-4911-3136-Gs.txt", Gs, ' ')