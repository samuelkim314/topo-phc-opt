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
chern =  parse(Int64, ARGS[3])

nin = 4
D = 3
res = 16
nbands = 10
nkpoints = 300
g = 12
random_init = true
objective = calccompletebandgap

# Create directories for saving input/output/log files
tmpdir = haskey(ENV, "TMPDIR")
sub_dir = "sg$sgnum-dim3-n$(@sprintf("%.2f", nin))/chern$chern/global-$id"
in_dir, out_dir, log_dir, in_dir_full, out_dir_full, out_dir_best_full = makedirs(sub_dir, tmpdir)

# Other parameters
cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)
# Has inversion symmetry or not
hasinv = SymOperation{3}("-x,-y,-z") in spacegroup(sgnum, D).operations

runtype = "all" # we don't have TE/TM in 3D; do all polarizations
gyro_type = "detfix"
epsin = nin^2
epsout = 1.0

if sgnum == 75
    if chern == 1
        topoindex = [1]
    elseif chern == 2
        topoindex = [2]
    elseif chern == 3
        topoindex = [3]
    end
elseif sgnum == 168
    if chern == 1
        topoindex = [5]
    elseif chern == 2
        topoindex = [4]
    elseif chern == 3
        topoindex = [3]
    end
else
    error("invalid sgnum")
end

# Random initialization
x0, nlatticeparams, nflatparams = randinit(sgnum, D, hasinv; cntr=cntr)

opt_id = 1
tracker = OptTracker()

function f(x::Vector, grad::Vector, opt_id::Integer)

    ff, pRs, flat_temp = params2flat(x, sgnum; 
        cntr=cntr, hasinv=hasinv, abclow=0.75, abchigh=1.25)
    Rs = conventionalize(pRs, cntr)

    # Symmetry calculations using MPB
    symcalcname = joinpath(sub_dir, "$opt_id-sym")
    band_topo = run_mpb_sym(symcalcname,
            (io;kwargs...)->prepare_mpbcalc_trsb!(io, sgnum, flat_temp, pRs, ff, epsin, epsout, runtype, 
                g, gyro_type;
                res=res, nbands=nbands, kwargs...),
            ()->checktopoindex(sgnum, symcalcname, topoindex=topoindex, dir=out_dir),
            false;
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

        # Note that mpbi requires both inversion symmetry and time-reversal
        # symmetry. We have broken the latter here, so we use regular mpb.
        bandgap, imax = run_mpb(shcalcname,
                        io->prepare_mpbcalc_trsb!(io, sgnum, flat_temp, pRs, ff, epsin, epsout, runtype,
                            g, gyro_type;
                            res=res, kvs=k_interp, id=id, nbands=nbands),
                        ()->objective(joinpath(sub_dir, calcname), band_topo, dir=out_dir);
                        hasinv=false,
                        dir=in_dir,
                        tmpdir=tmpdir)
        println(band_topo[imax])
    else
        bandgap = -2.0
    end

    return bandgap
end


# Perform optimization
tracker = globalopt_levelset(f, x0, nlatticeparams, nflatparams, hasinv, 5000,
    in_dir_full, out_dir_full, out_dir_best_full)