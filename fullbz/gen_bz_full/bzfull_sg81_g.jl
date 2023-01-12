using Crystalline
using MPBUtils
include("../../write_ctl.jl")

write_to_file = true

# load a previously generated structure
sgnum    = 81
D        = 3
subdir = "sg81-dim3-n4.00/ff/local-7-553-best"
loadstructname = joinpath(subdir, "1621")
# loadstructname = "thomas-candidates/dim3-sg147/1to6-14674-cladding-g0/local-2-7746-8732-cycle20-best/1844"
pRs, pflat, isoval, epsin, g, epsout, _ = lattice_from_mpbparams_trsb("output/$(loadstructname).sh") # Rs, flat, isoval, epsin, epsout, (irrfbz) kvecs 

# using PlotlyJS
# cntr = centering(sgnum, 3) # will be ‘C’ for base-centered
# Rs = conventionalize(pRs, cntr)
# kp = irrfbz_path(sgnum, Rs)
# println(kp)
# Gs = reciprocalbasis(Rs) # using Bravais to create the reciprocal basis
# cell = wignerseitz(Gs)   # construct associated Brillouin zone
# p = plot(latticize(cell, Gs))
# p2 = plot(latticize(kp, Gs))
# plotdata = append!(p.plot.data, p2.plot.data)
# PlotlyJS.plot(plotdata)
# pg = pointgroup(primitivize(spacegroup(sgnum, Val(3))))
# println([rotation(pgi) for pgi in pg])
# hasinv = SymOperation{3}("-x,-y,-z") in spacegroup(sgnum, D).operations
# println(spacegroup(sgnum, D).operations)


# computational parameters
nbands = 12
res = 32
runtype = "all"

# dense k-coverage of positive-kᵢ quadrant of the "plain" BZ
kvsearch = "bzsymreduced"
Nk = 40
if kvsearch == "bzfull"
    kᵢs = range(0.0, 1.0, length=Nk)
    kvs = [[kx, ky, kz] for kx in kᵢs, ky in kᵢs, kz in kᵢs]
elseif kvsearch == "bzsymreduced"
    kxs = range(0.0, 0.5, length=Nk)
    kys = range(0.0, 0.5, length=Nk)
    kzs = range(-0.5, 0.5, length=2*Nk+1)
    kvs = [[kx, ky, kz] for kx in kxs, ky in kys, kz in kzs]
end

# write to .sh file for mpb
if write_to_file
    write_dir = "output/"
    filename = loadstructname*"-bzsymreduced"
    kvecs_filename = subdir*"/$(kvsearch)_"*(contains(kvsearch, "bzsymreduced") ? "sg$(sgnum)_" : "")*"Nk$(Nk).dat"
    # kvecs_filename = "thomas-candidates/dim3-sg147/1to6-14674-cladding-g0/local-2-7746-8732-cycle20-best/$(kvsearch)_"*(contains(kvsearch, "bzsymreduced") ? "sg$(sgnum)_" : "")*"Nk$(Nk).dat"
    open(write_dir*filename*".sh", "w") do io
        prepare_mpbcalc_trsb!(io, sgnum, pflat, pRs, nothing, epsin, epsout, runtype, g;
                             res=res, kvecs=kvecs_filename, nbands=nbands, isoval=isoval)
    end
    println("MPB setup file written to:\n  ∙  ", write_dir*filename, ".sh\n")
    
    # write a k-vec file
    kvecs_path = write_dir*kvecs_filename
    if !isfile(kvecs_path)
        open(kvecs_path, "w") do io
            write(io, "(")
            foreach(kv->(write(io, "("); join(io, kv, " "); write(io, ") ")), kvs)
            write(io, ")")
        end
        println("New kvecs file written to:\n  ∙  ", kvecs_path)
    end
end
nothing