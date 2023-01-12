#=
Process NLopt results for optimizing bandgap.
=#
using MPBUtils
using DelimitedFiles
include("../bands.jl")
using Printf

nin = 4
sgnum = 147
sgnum_hinge = 147
D = 3
ntopobands = 6
g=0

prefix = "thomas-candidates/dim3-sg$sgnum/1to$ntopobands-14674-cladding-g$g"
# saveprefix = prefix
saveprefix = "local"

topopath = joinpath("thomas-candidates/dim3-sg$sgnum_hinge", 
        "1to$ntopobands-local-14674-best/2962")
gapbounds = getgap(topopath, ntopobands, dir="output")
gap0lower = gapbounds[1]
gap0upper = gapbounds[2]
println((gap0lower, gap0upper))
println(calcminbandgap(topopath, ntopobands))

topofullpath = "thomas-candidates/dim3-sg$sgnum_hinge/1to$ntopobands-local-14674-best/2962-bzsymreduced"
gapbounds = getgap(topofullpath, 4; dir="output")
gap0lower = gapbounds[1]
gap0upper = gapbounds[2]
println((gap0lower, gap0upper))
println(calcminbandgap(topofullpath, 4))
min, i = calcargminbandgap(topofullpath, ntopobands)

disp_data = readdlm(joinpath("output", topofullpath*"-dispersion.out"),',')
println(disp_data[i, :])

claddingpath = joinpath("thomas-candidates/dim3-sg$sgnum/",
        "1to$ntopobands-14674-cladding-g$g/local-2-7746-8732-cycle20-best/1844")
gapbounds = getgap(claddingpath, 7; dir="output")
gap0lower = gapbounds[1]
gap0upper = gapbounds[2]
println((gap0lower, gap0upper))

claddingpath = joinpath("thomas-candidates/dim3-sg$sgnum/",
        "1to$ntopobands-14674-cladding-g$g/local-2-7746-8732-cycle20-best/1844-bzsymreduced")
gapbounds = getgap(claddingpath, 7; dir="output")
gap0lower = gapbounds[1]
gap0upper = gapbounds[2]
println((gap0lower, gap0upper))