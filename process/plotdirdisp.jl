# Plot dispersion diagrams for all the files in a directory
# Useful for after running an optimization and sorting the best ones,
# e.g. using process_nlopt_gap.jl

using MPBUtils
using Printf
include("../write_ctl.jl")
include("../topology.jl")
include("../bands.jl")
using PlotlyJS

sgnum = 13
nin = 4

D = 3
mu1 = minconnectivity(sgnum)

trial_ids = [4, 7, 9]
global_best_ids = [579, 1242, 1916]

for trial in [9]
    # These are sub-folders in /input/ and /output/, respectively
    # Name is something like "sg68-n3.48-trial1"
    
    # localid = global_best_ids[trial]
    # trialid = trial_ids[trial]
    in_dir = "sg$sgnum-dim$D-n$(@sprintf("%.2f", nin))/opt-global-$trial-best"
    # in_dir = "sg$sgnum-dim$D-n$(@sprintf("%.2f", nin))/opt-local-$trialid-$localid-best"
    out_dir = in_dir
    in_dir_full = joinpath("input", in_dir)
    out_dir_full = joinpath("output", out_dir)
 
    fileslist = readdir(out_dir_full)
    # Sort out the files that match the format "$id.sh", e.g. "123.sh"
    fileslist = [filename for filename in fileslist if occursin(r"^[0-9]+\-dispersion.out$", filename)]

    println(fileslist)

    for filename in fileslist
        fullpath = joinpath(in_dir_full, filename)
        id = parse(Int, split(filename, "-")[1])
        calcname = "$id"

        completegap = calccompletebandgap(joinpath(out_dir, calcname), mu1)
        mingap = calcminbandgap(joinpath(out_dir, calcname), mu1)
        gap_bounds = getgap(joinpath(out_dir, calcname), mu1)

        # Extract lattice vectors to extract k-points
        in_path = joinpath(in_dir_full, calcname * ".sh")
        pRs, pflat, isoval, epsin, epsout, kvecs = lattice_from_mpbparams(in_path)
        cntr = centering(sgnum, 3) # will be ‘C’ for base-centered
        Rs = conventionalize(pRs, cntr)

        # println([gap_bounds[1], gap_bounds[2]])

        p = plotbands(joinpath(out_dir, calcname), Rs, sg=sgnum, dim=3,
            mu1=minconnectivity(sgnum),
            # title="$localid-$trialid-$id     Mingap: $(@sprintf("%.4f", mingap))    Compgap: $(@sprintf("%.4f", completegap))",
            # title="$trialid-$localid-$id     Bandgap: $(@sprintf("%.4f", completegap))",
            title="$trial-$id     Bandgap: $(@sprintf("%.4f", completegap))",
            # savename=joinpath(out_dir_full, "bands-nodal-$id.svg"), 
            nkpoints=300, ylim=[0.22, 0.37],
            # range_fill=[gap_bounds[1], gap_bounds[2]],
            width=600, height=250,
            # width=600, height=350,
            savename=joinpath(out_dir_full, "bands-zoom-$id.svg"),
            range_fill=[gap_bounds[1], gap_bounds[2]], 
            # range_fill=[0.315462299523985, 0.32873040920883245], fillcolor="rgba(166,51,255,0.2)"     # 4-579-3944
            # range_fill=[0.3066278376865743, 0.33454], fillcolor="rgba(166,51,255,0.2)"  # 7-1242-3918
            # range_fill=[0.3222818052499316, 0.3648167060788395], fillcolor="rgba(166,51,255,0.2)"  # 9-1916-4417
            )
        # PlotlyJS.savefig(p, joinpath(out_dir_full, "bands-$id.png"))
        # display(p)
        # return p
        # println(p.plot.layout)
        
    end
end