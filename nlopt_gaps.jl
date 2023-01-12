# Optimize any gap with non-trivial bands below.
using Printf
include("write_ctl.jl")
include("bands.jl")
include("process_symeigs.jl")
include("topology.jl")
include("params.jl")
include("opt_utils.jl")
include("random.jl")
using Crystalline
using LinearAlgebra

id = parse(Int64, ARGS[1])
sgnum = parse(Int64, ARGS[2])
# 27, 37, 81, 82, 103 and 184 can achieve Weyl points with matched frequencies.
# 3, 75, 77, 168, 171, 172 are all bad.

nin = 4     # Material refractive index
D = 3       # Dimension (2D or 3D)
res = 16    # Simulation resolution
nbands = 8  # Number of bands to simulate. Will check up to nbands-1 for topology.
nkpoints = 300  # Number of k-points. ~100 is probably sufficient.
objective = calccompletebandgap

# Create directories for saving input/output/log files
tmpdir = haskey(ENV, "TMPDIR")
sub_dir = "sg$sgnum-dim3-n$(@sprintf("%.2f", nin))/gaps/global-$id"
in_dir, out_dir, log_dir, in_dir_full, out_dir_full, out_dir_best_full = makedirs(sub_dir, tmpdir)

# Other parameters
cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)
# Has inversion symmetry or not
hasinv = SymOperation{3}("-x,-y,-z") in spacegroup(sgnum, D).operations

epsin = nin^2
epsout = 1.0
runtype = "all" # we don't have TE/TM in 3D; do all polarizations

x0, nlatticeparams, nflatparams = randinit(sgnum, D, hasinv; cntr=cntr)


function f(x::Vector, grad::Vector, opt_id::Integer)
    ff, pRs, flat_temp = params2flat(x, sgnum; 
        cntr=cntr, hasinv=hasinv, abclow=0.75, abchigh=1.25)
    Rs = conventionalize(pRs, cntr)

    # Symmetry calculations using MPB
    symcalcname = joinpath(sub_dir, "$opt_id-sym")
    band_topo = run_mpb_sym(symcalcname,
            (io;kwargs...)->prepare_mpbcalc!(io, sgnum, flat_temp, pRs, ff, epsin, epsout, runtype;
                res=res, nbands=nbands, kwargs...),
            ()->checknontriviallowest(sgnum, symcalcname, dir=out_dir),
            true;
            indir=in_dir,
            tmpdir=tmpdir
        )

    # Run MPB simulation for full dispersion
    if !isempty(band_topo)
        # k-points
        kp = irrfbz_path(sgnum, Rs)
        k_interp = interpolate(kp, nkpoints)

        calcname = "$opt_id"
        shcalcname = joinpath(sub_dir, calcname)
        bands_list = collect(keys(band_topo))

        bandgap, imax = run_mpb(shcalcname,
                        io->prepare_mpbcalc!(io, sgnum, flat_temp, pRs, ff, epsin, epsout, runtype;
                                res=res, kvs=k_interp, id=id, nbands=nbands),
                        ()->objective(joinpath(sub_dir, calcname), bands_list, dir=out_dir);
                        hasinv=hasinv,
                        dir=in_dir,
                        tmpdir=tmpdir)
        println(bands_list[imax])
        println(band_topo[bands_list[imax]])
    else
        bandgap = -2.0
    end

    return bandgap

end

# Perform optimization
tracker = globalopt(f, x0, nlatticeparams, nflatparams, hasinv, 5000,
    in_dir_full, out_dir_full, out_dir_best_full)