# Plot dispersion diagrams for all the files in a directory
# Useful for after running an optimization and sorting the best ones,
# e.g. using process_nlopt_gap.jl
using MPBUtils
using Printf
include("../src/write_ctl.jl")
using PyPlot # This is required for mesh_3d_levelsetlattice
using Crystalline
using GLMakie
include("plot_makie.jl")

sgnum = 75
nin = 4
D = 3

# global_best_ids = [839, 924, 446, 955, 461]
# global_best_ids = [759, 775, 894, 303, 643]
# local_best_ids = [336, 292, 273, 389, 210]

trial_ids = [3]
global_best_ids = [3186]
# trial_ids = [1, 2, 4, 8, 10]    # SG 82
# global_best_ids = [4286, 4911, 3886, 4882, 4548]

cntr = centering(sgnum, D)# will be ‘C’ for base-centered

for trial in [1]
    # These are sub-folders in /input/ and /output/, respectively
    # Name is something like "sg68-n3.48-trial1"
    
    localid = global_best_ids[trial]
    trialid = trial_ids[trial]
    # in_dir = "sg$sgnum-dim$D-n$(@sprintf("%.2f", nin))/opt-global-$trial-best"
    in_dir = "sg$sgnum-dim$D-n$(@sprintf("%.2f", nin))/chern2-ff-direct/local-$trialid-$localid-best"
    out_dir = in_dir
    in_dir_full = joinpath("output", in_dir)
    out_dir_full = joinpath("output", out_dir)
 
    fileslist = readdir(out_dir_full)
    # Sort out the files that match the format "$id.sh", e.g. "123.sh"
    fileslist = [filename for filename in fileslist if occursin(r"^[0-9]+\-dispersion.out$", filename)]

    println(fileslist)

    for filename in fileslist
        fullpath = joinpath(in_dir_full, filename)
        id = parse(Int, split(filename, "-")[1])
        calcname = "$id"

        # Extract lattice vectors to extract k-points
        in_path = joinpath(in_dir_full, calcname * ".sh")
        pRs, pflat, isoval, epsin, epsout, kvecs = lattice_from_mpbparams(in_path)
        Rs = conventionalize(pRs, cntr)
        flat = conventionalize(pflat, cntr)

        # f = figure()
        # Crystalline.plot(flat, Rs; isoval=isoval, fig=f)
        # PyPlot.savefig(joinpath(out_dir_full,"pc-$id.png"))
        # close(f)
        # return f

        # f = plot_pc_rotated(flat, Rs, isoval, 64)
        f = plot_pc(flat, Rs, isoval, 64)
        GLMakie.save(joinpath(out_dir_full, "pc-$id.png"), f)
        # display(f)
    end
end