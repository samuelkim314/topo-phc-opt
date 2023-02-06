# Global bandgap optimization of Gamma-enforced (filling-enforced) gap
using Printf
include("src/write_ctl.jl")
include("src/bands.jl")
include("src/process_symeigs.jl")
include("src/topology.jl")
include("src/params.jl")
include("src/opt_utils.jl")
include("src/random.jl")
using Crystalline
using LinearAlgebra

id = parse(Int64, ARGS[1])
sgnum = parse(Int64, ARGS[2])

nin = 4
D = 3
res = 16
nbands = 8
nkpoints = 300
random_init = true
objective = calccompletebandgap

# Create directories for saving input/output/log files
tmpdir = haskey(ENV, "TMPDIR")
sub_dir = "sg$sgnum-dim3-n$(@sprintf("%.2f", nin))/filling/global-$id"
in_dir, out_dir, log_dir, in_dir_full, out_dir_full, out_dir_best_full = makedirs(sub_dir, tmpdir)


# Other parameters
mu1 = minconnectivity(sgnum)
cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)
# Has inversion symmetry or not
hasinv = SymOperation{3}("-x,-y,-z") in spacegroup(sgnum, D).operations
epsin = nin^2
epsout = 1.0
runtype = "all" # we don't have TE/TM in 3D; do all polarizations

x0, nlatticeparams, nflatparams = randinit(sgnum, D, hasinv; cntr=cntr)


function f(x::Vector, grad::Vector, opt_id::Integer)
    # Extract lattice/geometry parameters
    ff, pRs, flat_temp = params2flat(x, sgnum; 
        cntr=cntr, hasinv=hasinv, abclow=0.75, abchigh=1.25)
    Rs = conventionalize(pRs, cntr)

    # Symmetry calculations using MPB
    symcalcname = joinpath(sub_dir, "$opt_id-sym")
    band_topo = run_mpb_sym(symcalcname,
            (io;kwargs...)->prepare_mpbcalc!(io, sgnum, flat_temp, pRs, ff, epsin, epsout, runtype;
                res=res, nbands=nbands, kwargs...),
            ()->calcbandtopo(sgnum, symcalcname, dir=out_dir),
            true;
            indir=in_dir,
            tmpdir=tmpdir
        )

    # Run MPB simulation for full dispersion
    # Only look at bottom bands.
    if haskey(band_topo, 1:mu1) && band_topo[1:mu1] != TRIVIAL
        # k-points
        kp = irrfbz_path(sgnum, Rs)
        k_interp = interpolate(kp, nkpoints)

        calcname = "$opt_id"
        shcalcname = joinpath(sub_dir, calcname)

        bandgap = run_mpb(shcalcname,
                        io->prepare_mpbcalc!(io, sgnum, flat_temp, pRs, ff, epsin, epsout, runtype;
                                res=res, kvs=k_interp, id=id, nbands=nbands),
                        ()->objective(joinpath(sub_dir, calcname), mu1, dir=out_dir);
                        hasinv=hasinv,
                        dir=in_dir,
                        tmpdir=tmpdir)
    else
        if objective == calcminbandgap
            bandgap = 0.0
        else
            bandgap = -2.0
        end
    end

    return bandgap
end

# Perform optimization
tracker = globalopt(f, x0, nlatticeparams, nflatparams, hasinv, 5000,
    in_dir_full, out_dir_full, out_dir_best_full)