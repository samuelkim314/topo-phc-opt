using NLopt
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
sgnum = parse(Int64, ARGS[2]) # 81 or 147
chern =  parse(Int64, ARGS[3])

nin = 4
D = 3
res = 16
nbands = 10
nkpoints = 300
g = 12
random_init = false
objective = calccompletebandgap

# Optimization parameters
iters_per_cycle = 500
ncycles = 10

# Initial conditions
if sgnum == 75
    if chern == 1
        topoindex = [1]
        global_best_trials = [1, 3, 4, 6, 7, 8]
        global_best_ids = [4937, 43, 3780, 4706, 4964, 4919]
        global_best_trials = [12, 14, 17, 19]
        global_best_ids = [3762, 3898, 2300, 2419]
    elseif chern == 2
        topoindex = [2]
        global_best_trials = [1, 3, 4, 7, 8, 9]
        global_best_ids = [497, 3186, 224, 4754, 4973, 4742]
    elseif chern == 3
        topoindex = [3]
        global_best_trials = [5]
        global_best_ids = [4958]
    end
elseif sgnum == 168
    if chern == 1
        topoindex = [5]
        global_best_trials = [1, 3, 5, 7, 9]
        global_best_ids = [4927, 4912, 4611, 3329, 4642]
    elseif chern == 2
        topoindex = [4]
        global_best_trials = [2, 4, 6, 10]
        global_best_ids = [1455, 390, 4634, 2790]
    end
else
    error("invalid sgnum")
end
trialid = global_best_trials[id]
localid = global_best_ids[id]


# Create directories for saving input/output/log files
tmpdir = haskey(ENV, "TMPDIR")
sub_dir = "sg$sgnum-dim3-n$(@sprintf("%.2f", nin))/chern$chern/local-$trialid-$localid"
in_dir, out_dir, log_dir, in_dir_full, out_dir_full, out_dir_best_full = makedirs(sub_dir, tmpdir)

# Other parameters
cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)
# Has inversion symmetry or not
hasinv = SymOperation{3}("-x,-y,-z") in spacegroup(sgnum, D).operations


runtype = "all" # we don't have TE/TM in 3D; do all polarizations
gyro_type = "detfix"
epsin = nin^2
epsout = 1.0

# Getting the initial condition
init_path = "output/sg$sgnum-dim3-n$(@sprintf("%.2f", nin))/chern$chern/global-$trialid-best/$localid.sh"
x0, g, nlatticeparams, nflatparams = sh_init_trsb(init_path, sgnum, hasinv; cntr=cntr, abclow=0.75, abchigh=1.25)


function f(x::Vector, grad::Vector, opt_id::Integer)
    ff, pRs, flat_temp = params2flat(x, sgnum; cntr=cntr, hasinv=hasinv, abclow=0.75, abchigh=1.25)
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
    if !isnothing(band_topo) && !isempty(band_topo)
        # k-points
        kp = irrfbz_path(sgnum, Rs)
        k_interp = interpolate(kp, nkpoints)
        calcname = "$opt_id"
        shcalcname = joinpath(sub_dir, calcname)

        # Run MPB and get the objective
        bandgap, _ = run_mpb(shcalcname,
                        io->prepare_mpbcalc_trsb!(io, sgnum, flat_temp, pRs, ff, epsin, epsout, runtype,
                            g, gyro_type;
                            res=res, kvs=k_interp, id=id, nbands=nbands),
                        ()->objective(joinpath(sub_dir, calcname), band_topo; dir=out_dir);
                        hasinv=false,
                        dir=in_dir,
                        tmpdir=tmpdir)
    else
        bandgap = -2.0
    end

    return bandgap

end

# Perform optimization!
localopt_levelset(f, x0, nlatticeparams, nflatparams, hasinv, 
    ncycles, iters_per_cycle, in_dir_full, out_dir_full, out_dir_best_full)
