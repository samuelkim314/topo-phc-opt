# Optimize any gap with non-trivial bands below.
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
objective = calccompletebandgap

# Optimization parameters
iters_per_cycle = 500
ncycles = 10

if sgnum==27
    trial_ids = [1, 2, 6, 7, 8, 10]
    global_best_ids = [4850, 4827, 4724, 484, 434, 482]
    trial_ids = [1, 2, 3, 4, 5]
    global_best_ids = [6, 101, 41, 11, 30]
elseif sgnum==37
    # thomas opt-global-ff-direct
    trial_ids = [1, 2, 6, 9]
    global_best_ids = [4634, 517, 904, 4279]
elseif sgnum==75
    trial_ids = [1, 3, 9]
    global_best_ids = [4483, 4092, 455]
elseif sgnum==81
    trial_ids = [1, 7, 9, 10]
    global_best_ids = [4210, 349, 3098, 4822]
elseif sgnum==82
    trial_ids = [1, 2, 4, 8, 10]
    global_best_ids = [4286, 4911, 3886, 4882, 4548]
elseif sgnum==184
    trial_ids = [4, 6, 7, 8, 9]
    global_best_ids = [4951, 4986, 3315, 807, 4944]
elseif sgnum==103
    trial_ids = [5, 6, 8, 9, 10]
    global_best_ids = [3097, 933, 913, 4625, 3527]
end
trial = trial_ids[id]
globalid = global_best_ids[id]


# Create directories for saving input/output/log files
tmpdir = haskey(ENV, "TMPDIR")
sub_dir = "sg$sgnum-dim3-n$(@sprintf("%.2f", nin))/gaps/local-$trial-$globalid"
in_dir, out_dir, log_dir, in_dir_full, out_dir_full, out_dir_best_full = makedirs(sub_dir, tmpdir)

# Other parameters
cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)
# Has inversion symmetry or not
hasinv = SymOperation{3}("-x,-y,-z") in spacegroup(sgnum, D).operations

runtype = "all" # we don't have TE/TM in 3D; do all polarizations
epsin = nin^2
epsout = 1.0

# Getting the initial condition
init_path = "output/sg$sgnum-dim3-n$(@sprintf("%.2f", nin))/gaps/global-$trial-best/$globalid.sh"
x0, nlatticeparams, nflatparams = sh_init(init_path, sgnum, hasinv; cntr=cntr, abclow=0.75, abchigh=1.25)

function f(x::Vector, grad::Vector)
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
                        tmpdir=true)
        println(bands_list[imax])
        println(band_topo[bands_list[imax]])
    else
        bandgap = -2.0
    end

    return bandgap

end

# Perform optimization!
localopt_levelset(f, x0, nlatticeparams, nflatparams, hasinv, ncycles, iters_per_cycle,
    in_dir_full, out_dir_full, out_dir_best_full)
