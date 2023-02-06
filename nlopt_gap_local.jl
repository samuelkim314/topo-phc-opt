# Local bandgap optimization of Gamma-enforced (filling-enforced) gap
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
sgnum = parse(Int64, ARGS[2])   # 13, 48, 49, 50, 68, 86

nin = 4
D = 3
res = 16
nbands = 8
nkpoints = 300
random_init = false
objective = calccompletebandgap

# Optimization parameters
iters_per_cycle = 500
ncycles = 10

if sgnum==68
    # trial_ids = [2, 2, 3, 3, 3, 5, 6, 6, 7]
    # global_best_ids = [342, 1167, 1455, 1883, 1143, 1248, 1617, 844, 1724]
    
    # opt-global-wideisoval
    trial_ids = [1, 3, 3, 4, 5, 6, 7, 8, 20]
    global_best_ids = [9681, 6089, 8134, 7599, 9813, 9801, 7083, 8930, 6350]
elseif sgnum==13
    trial_ids = [1, 2, 2, 3, 4, 4, 6, 7, 7, 9]
    global_best_ids = [1560, 1582, 1736, 445, 579, 1587, 1804, 1242, 1572, 1916]
    trial_ids = [3, 5]
    global_best_ids = [989, 4997]
elseif sgnum==48
    trial_ids = [1, 2, 3, 4, 5]
    global_best_ids = [1018, 4297, 12, 4915, 488]
elseif sgnum==49
    trial_ids = [2, 2, 2, 3, 3, 8, 10]
    global_best_ids = [805, 677, 107, 1860, 860, 933, 1151]
    
elseif sgnum==50
    trial_ids = [1, 2, 5, 5, 6, 7, 7, 8]
    global_best_ids = [1777, 683, 1132, 1055, 1785, 1251, 1558, 534]

elseif sgnum==86
    trial_ids = [2, 2, 4, 4, 8, 8, 8, 8, 8, 9]
    global_best_ids = [1541, 1466, 1645, 1528, 1307, 1545, 1845, 1457, 1169, 1532]
end
trial = trial_ids[id]
globalid = global_best_ids[id]


# Create directories for saving input/output/log files
tmpdir = haskey(ENV, "TMPDIR")
sub_dir = "sg$sgnum-dim3-n$(@sprintf("%.2f", nin))/filling/local2-$trial-$globalid"
in_dir, out_dir, log_dir, in_dir_full, out_dir_full, out_dir_best_full = makedirs(sub_dir, tmpdir)


# Other parameters
mu1 = minconnectivity(sgnum)
cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)
# Has inversion symmetry or not
hasinv = SymOperation{3}("-x,-y,-z") in spacegroup(sgnum, D).operations

runtype = "all" # we don't have TE/TM in 3D; do all polarizations
epsin = nin^2
epsout = 1.0

# Getting the initial condition
init_path = "output/sg$sgnum-dim3-n$(@sprintf("%.2f", nin))/filling/global-$trial-best/$globalid.sh"
x0, nlatticeparams, nflatparams = sh_init(init_path, sgnum, hasinv; cntr=cntr, abclow=0.75, abchigh=1.25)


function f(x::Vector, grad::Vector, opt_id::Integer)

    ff, pRs, flat_temp = params2flat(x, sgnum; cntr=cntr, hasinv=hasinv, abclow=0.75, abchigh=1.25)
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

    # Only look at bottom bands.
    # Run MPB simulation for full dispersion
    if haskey(band_topo, 1:mu1) && band_topo[1:mu1] != TRIVIAL
        # k-points
        kp = irrfbz_path(sgnum, Rs)
        k_interp = interpolate(kp, nkpoints)
        calcname = "$opt_id"
        shcalcname = joinpath(sub_dir, calcname)

        # Run MPB and get the objective
        bandgap = run_mpb(shcalcname,
                        io->prepare_mpbcalc!(io, sgnum, flat_temp, pRs, ff, epsin, epsout, runtype;
                                res=res, kvs=k_interp, id=id, nbands=nbands),
                        ()->objective(joinpath(sub_dir, calcname), 2; dir=out_dir);
                        hasinv=hasinv,
                        dir=in_dir,
                        tmpdir=tmpdir)
    else
        bandgap = -2.0
    end

    return bandgap

end

# Perform optimization!
localopt(f, x0, nlatticeparams, nflatparams, hasinv, ncycles, iters_per_cycle,
    in_dir_full, out_dir_full, out_dir_best_full)