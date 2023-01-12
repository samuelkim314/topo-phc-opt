using Crystalline
using MPBUtils
include("../../write_ctl.jl")

write_to_file = true

# load a previously generated structure
sgnum    = 147
D        = 3
loadstructname = "thomas-candidates/dim3-sg147/1to6-local-14674-best/2962"
# loadstructname = "thomas-candidates/dim3-sg147/1to6-14674-cladding-g0/local-2-7746-8732-cycle20-best/1844"
subdir = "sg147-dim3-n4.00/ff-nbands18/local-2-2390-best"
loadstructname = joinpath(subdir, "4321")
pRs, pflat, isoval, epsin, g, epsout, _ = lattice_from_mpbparams_trsb("output/$(loadstructname).sh") # Rs, flat, isoval, epsin, epsout, (irrfbz) kvecs 

# computational parameters
nbands = 16
res = 32
runtype = "all"

# dense k-coverage of positive-kᵢ quadrant of the "plain" BZ
kvsearch = "bzsymreduced"
Nk = 40
if kvsearch == "bzfull"
    kᵢs = range(0.0, 1.0, length=Nk)
    kvs = [[kx, ky, kz] for kx in kᵢs, ky in kᵢs, kz in kᵢs]
elseif kvsearch == "bzsymreduced"
    kxs = range(0, 0.5, length=Nk)
    kys = range(0.0, 1/3,  length=Nk)
    kzs = range(-0.5, 0.5,  length=2*Nk+1)
    kvs = [[kx, ky, kz] for kx in kxs, ky in kys, kz in kzs]
end

# write to .sh file for mpb
if write_to_file
    # write_dir = "input/"
    write_dir = "output/"
    filename = loadstructname*"-bzsymreduced"
    kvecs_filename = subdir*"/$(kvsearch)_"*(contains(kvsearch, "bzsymreduced") ? "sg$(sgnum)_" : "")*"Nk$(Nk).dat"
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