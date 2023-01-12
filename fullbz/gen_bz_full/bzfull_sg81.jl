using Crystalline
using MPBUtils

write_to_file = true

# load a previously generated structure
sgnum    = 81
D        = 3
loadstructname = "sg$sgnum-dim3-n4.00/gaps-ff-direct/local-9-3098-best/4990"
loadstructname = "sg$sgnum-dim3-n4.00/gaps-ff-direct/local-10-4822-best/4358"
loadstructname = "sg$sgnum-dim3-n4.00/gaps-ff-direct/local-9-3098-best/4384"
loadstructname = "sg$sgnum-dim3-n4.00/gaps-ff-direct/local-7-349-best/2840"
pRs, pflat, isoval, epsin, epsout, _ = lattice_from_mpbparams("output/$(loadstructname).sh") # Rs, flat, isoval, epsin, epsout, (irrfbz) kvecs 

# using PlotlyJS
# cntr = centering(sgnum, 3) # will be ‘C’ for base-centered
# Rs = conventionalize(pRs, cntr)
# kp = irrfbz_path(sgnum, Rs)
# println(kp.points)
# Gs = reciprocalbasis(Rs) # using Bravais to create the reciprocal basis
# cell = wignerseitz(Gs)   # construct associated Brillouin zone
# p = PlotlyJS.plot(cell)
# p2 = PlotlyJS.plot(kp)
# plotdata = append!(p.plot.data, p2.plot.data)
# p = PlotlyJS.plot(plotdata)
# display(p)
# p = PlotlyJS.plot(latticize(kp, Rs))
# display(p)
# asdf

# computational parameters
nbands = 8
res = 32
runtype = "all"

# dense k-coverage of positive-kᵢ quadrant of the "plain" BZ
kvsearch = "bzsymreduced"
Nk = 80
if kvsearch == "bzfull"
    kᵢs = range(0.0, 1.0, length=Nk)
    kvs = [[kx, ky, kz] for kx in kᵢs, ky in kᵢs, kz in kᵢs]
elseif kvsearch == "bzsymreduced"
    kxs = range(0.0, 0.5, length=Nk)
    kys = range(0.0, 0.5, length=Nk)
    kzs = range(0.0, 0.5, length=Nk)
    kvs = [[kx, ky, kz] for kx in kxs, ky in kys, kz in kzs]
elseif kvsearch == "bzsymreduced_nodalzoom116"
    kxs = range(0.4, 0.5,       length=Nk÷5)
    kys = range(-0.065, 0.135,  length=(Nk*2)÷5)
    kzs = range(0.0, 0.5,       length=Nk)
    kvs = [[kx, ky, kz] for kx in kxs, ky in kys, kz in kzs]
end

# write to .sh file for mpb
if write_to_file
    write_dir = "output/"
    filename = loadstructname*"-bzsymreduced"
    kvecs_filename = "sg$sgnum-dim3-n4.00/gaps-ff-direct/local-9-3098-best/$(kvsearch)_"*(contains(kvsearch, "bzsymreduced") ? "sg$(sgnum)_" : "")*"Nk$(Nk).dat"
    open(write_dir*filename*".sh", "w") do io
        prepare_mpbcalc!(io, sgnum, pflat, pRs, nothing, epsin, epsout, runtype;
                             res=res, kvecs=kvecs_filename, nbands=nbands, isoval=isoval)
    end
    println("MPB setup file written to:\n  ∙  ", write_dir*filename, ".sh\n")
    
    # write a k-vec file
    kvecs_path = write_dir*kvecs_filename
    open(kvecs_path, "w") do io
        write(io, "(")
        foreach(kv->(write(io, "("); join(io, kv, " "); write(io, ") ")), kvs)
        write(io, ")")
    end
    println("New kvecs file written to:\n  ∙  ", kvecs_path)

end
nothing