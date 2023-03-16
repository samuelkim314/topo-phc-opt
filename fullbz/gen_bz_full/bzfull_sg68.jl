using Crystalline, PyPlot
using MPBUtils

write_to_file = true

# load a previously generated structure
sgnum    = 68
D        = 3
resload  = 32
calcidx  = 116
loadstructname = "dim$(D)-sg$(sgnum)-irrfbz_$(calcidx)_Nk100-res$(resload)"
pRs, pflat, isoval, epsin, epsout, _ = lattice_from_mpbparams((@__DIR__)*"/../../input/$(loadstructname).sh") # Rs, flat, isoval, epsin, epsout, (irrfbz) kvecs 

# computational parameters
nbands = 4
res = 32
runtype = "all"
#plot(pflat, pRs; isoval=isoval)

# dense k-coverage of positive-kᵢ quadrant of the "plain" BZ
kvsearch = "bzsymreduced_nodalzoom116"
Nk = 110
if kvsearch == "bzfull"
    kᵢs = range(0.0, 1.0, length=Nk)
    kvs = [[kx, ky, kz] for kx in kᵢs, ky in kᵢs, kz in kᵢs]
elseif kvsearch == "bzsymreduced"
    kxs = range(-0.5, 0.5, length=2Nk)
    kys = range(0.0, 0.5,  length=Nk)
    kzs = range(0.0, 0.5,  length=Nk)
    kvs = [[kx, ky, kz] for kx in kxs, ky in kys, kz in kzs]
end

# write to .sh file for mpb
id = "$(kvsearch)_$(calcidx)_Nk$(Nk)"
if write_to_file
    write_dir = (@__DIR__)*"/../../input/"
    kvecs_filename = "$(kvsearch)_"*(contains(kvsearch, "bzsymreduced") ? "sg$(sgnum)_" : "")*"Nk$(Nk).dat"
    filename = mpb_calcname(D, sgnum, id, res, runtype)
    open(write_dir*filename*".sh", "w") do io
        prepare_mpbcalc!(io, sgnum, pflat, pRs, nothing, epsin, epsout, runtype;
                             res=res, kvecs=kvecs_filename, id=id, nbands=nbands, isoval=isoval)
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