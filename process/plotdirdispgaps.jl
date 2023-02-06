# Plot dispersion diagrams for all the files in a directory
# Useful for after running an optimization and sorting the best ones,
# e.g. using process_nlopt_gap.jl

using MPBUtils
using Printf
include("../src/write_ctl.jl")
include("../src/topology.jl")
include("../src/bands.jl")
include("../src/process_symeigs.jl")
# using PlotlyJS

sgnum = 82
nin = 4
D = 3

cntr = centering(sgnum, 3) # will be ‘C’ for base-centered

trial_ids = [1, 2, 4, 8, 10]    # SG 82
global_best_ids = [4286, 4911, 3886, 4882, 4548]
# trial_ids = [1, 2, 6, 7, 8, 10]     # SG 27
# global_best_ids = [4850, 4827, 4724, 484, 434, 482]
# trial_ids = [6]
# global_best_ids = [904]
# trial_ids = [9, 10]
# global_best_ids = [3098, 4822]
# trial_ids = [7, 9]
# global_best_ids = [349, 3098]
# trial_ids = [4]
# global_best_ids = [3886]

for trial in [1]
    # These are sub-folders in /input/ and /output/, respectively
    # Name is something like "sg68-n3.48-trial1"
    
    localid = global_best_ids[trial]
    trialid = trial_ids[trial]
    # in_dir = "sg$sgnum-dim$D-n$(@sprintf("%.2f", nin))/opt-global-$trial-best"
    in_dir = "sg$sgnum-dim$D-n$(@sprintf("%.2f", nin))/ff-direct/rand-local-$trialid-$localid-best"
    out_dir = in_dir
    in_dir_full = joinpath("output", in_dir)
    out_dir_full = joinpath("output", out_dir)
 
    fileslist = readdir(out_dir_full)
    # Sort out the files that match the format "$id.sh", e.g. "123.sh"
    fileslist = [filename for filename in fileslist if occursin(r"^[0-9]+\-dispersion.out$", filename)]

    println(fileslist)

    # for filename in fileslist
    for filename in ["4997-16bands-dispersion.out"]
        fullpath = joinpath(in_dir_full, filename)
        id = parse(Int, split(filename, "-")[1])
        calcname = "$id"

        band_topo = checknontriviallowest(sgnum, joinpath(out_dir, calcname*"-sym"), dir="output")
        bands_list = collect(keys(band_topo))

        completegap, imax = calccompletebandgap(joinpath(out_dir, calcname), bands_list, dir="output")
        mingap, _ = calcminbandgap(joinpath(out_dir, calcname), bands_list, dir="output")
        topoindex = band_topo[bands_list[imax]]
        bandsbelow = bands_list[imax][end]
        println(topoindex)

        gap_bounds = getgap(joinpath(out_dir, calcname), bandsbelow, dir="output")

        # Extract lattice vectors to extract k-points
        in_path = joinpath(out_dir_full, calcname * ".sh")
        pRs, pflat, isoval, epsin, epsout, kvecs = lattice_from_mpbparams(in_path)
        Rs = conventionalize(pRs, cntr)

        p = plotbands(joinpath(out_dir, calcname*"-16bands"), Rs, sg=sgnum, dim=3,
            mu1=bandsbelow,
            title="$sgnum-$trialid-$localid-$id    Compgap: $(@sprintf("%.4f", completegap))"*
                "    Topological index: $topoindex    Bands below: $bandsbelow",
            savename=joinpath(out_dir_full, "bands-zoom-$id.svg"), nkpoints=300,
            ylim=[0.6, 0.73], 
            range_fill=[gap_bounds[1], gap_bounds[2]], 
            # hlinesy=[0.53, 0.68],
            width=600, height=250,
            hlinesy=[0.6633]
            )

        display(p)
        # return p
       
    end
end