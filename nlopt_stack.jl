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

epsin1 = 11.68  # silicon
epsin2 = 3.9    # silica
D = 3
res = 16
nbands = 8
nkpoints = 300
random_init = true
objective = calccompletebandgap

# Create directories for saving input/output/log files
tmpdir = haskey(ENV, "TMPDIR")
sub_dir = "sg$sgnum-dim2/filling/global-$id"
in_dir, out_dir, log_dir, in_dir_full, out_dir_full, out_dir_best_full = makedirs(sub_dir, tmpdir)


# Other parameters
mu1 = 2     # Band index above which to maximize bandgap
cntr = centering(sgnum, D-1) # centering symbol (e.g, 'F' for face-centered, etc.)
# Has inversion symmetry or not
hasinv = SymOperation{3}("-x,-y,-z") in spacegroup(sgnum, D-1).operations
epsout = 1.0
runtype = "all" # we don't have TE/TM in 3D; do all polarizations

ngeomparams = 2

nlengthparams, nangleparams = nbasisparams(sgnum, Val(D-1))     # Wallpaper group
nlatticeparams = nlengthparams + nangleparams + 1   # Add one for c
x0 = rand(ngeomparams + nlatticeparams)


function f(x::Vector, grad::Vector, opt_id::Integer)
    # Extract lattice/geometry parameters
    basisparams = x[1:nlatticeparams-1]
    c = x[nlatticeparams]
    geomparams = x[nlatticeparams+1:end]
    Rs = params2basis(basisparams, sgnum, Val(D-1); abclow=0.75, abchigh=1.25)
    pRs = primitivize(Rs, cntr)
    Rs3D = extrudebasis(Rs, c)
    pRs3D = primitivize(Rs3D, cntr)

    # TODO: unnormalize the geometry parameters

    # Symmetry calculations using MPB
    symcalcname = joinpath(sub_dir, "$opt_id-sym")
    band_topo = run_mpb_sym(symcalcname,
            (io;kwargs...)->prepare_mpbcalc_geom!(io, sgnum, pRs3D, epsin1, epsin2, epsout, 
                geomparams[1], geomparams[2], runtype;
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
                        io->prepare_mpbcalc_geom!(io, sgnum, pRs3D, epsin1, epsin2, epsout, 
                            geomparams[1], geomparams[2], runtype;
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
tracker = globalopt(f, x0, 5000,
    in_dir_full, out_dir_full, out_dir_best_full)