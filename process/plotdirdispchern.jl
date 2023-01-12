# Plot dispersion diagrams for all the files in a directory
# Useful for after running an optimization and sorting the best ones,
# e.g. using process_nlopt_gap.jl

using MPBUtils
using Printf
include("../write_ctl.jl")
include("../topology.jl")
include("../bands.jl")
include("../process_symeigs.jl")

# sgnum = parse(Int64, ARGS[1])
sgnum = 168     # 75 or 168
chern = 1
nin = 4
D = 3

prefix = "local"

if sgnum == 75
    if chern == 1
        topoindex = [1]
        global_best_trials = [1, 3, 4, 6, 7, 8]
        global_best_ids = [4937, 43, 3780, 4706, 4964, 4919]
        
    elseif chern == 2
        topoindex = [2]
        global_best_trials = [1, 3, 4, 7, 8, 9]
        global_best_ids = [497, 3186, 224, 4754, 4973, 4742]
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

cntr = centering(sgnum, 3) # will be ‘C’ for base-centered

for trial in [3]
    # These are sub-folders in /input/ and /output/, respectively
    # Name is something like "sg68-n3.48-trial1"
    
    trialid = global_best_trials[trial]
    localid = global_best_ids[trial]
    in_dir = "sg$sgnum-dim$D-n$(@sprintf("%.2f", nin))/chern$chern-ff-direct/$prefix-$trialid-$localid-best"
    out_dir = in_dir
    in_dir_full = joinpath("output", in_dir)
    out_dir_full = joinpath("output", out_dir)
 
    fileslist = readdir(out_dir_full)
    # Sort out the files that match the format "$id.sh", e.g. "123.sh"
    fileslist = [filename for filename in fileslist if occursin(r"^[0-9]+\-dispersion.out$", filename)]

    println(fileslist)

    for filename in fileslist[1:end]
        fullpath = joinpath(in_dir_full, filename)
        id = parse(Int, split(filename, "-")[1])
        calcname = "$id"

        band_topo = checktopoindex(sgnum, joinpath(in_dir, calcname*"-sym"), topoindex=topoindex, dir="output")
        println(band_topo)
        completegap, imax = calccompletebandgap(joinpath(out_dir, calcname), band_topo)
        bandsbelow = band_topo[imax][end]
        mingap, imax = calcminbandgap(joinpath(out_dir, calcname), band_topo)

        gap_bounds = getgap(joinpath(out_dir, calcname), bandsbelow, dir="output")
        
        # Extract lattice vectors to extract k-points
        in_path = joinpath(in_dir_full, calcname * ".sh")
        pRs, pflat, isoval, epsin, epsout, kvecs = lattice_from_mpbparams(in_path)
        Rs = conventionalize(pRs, cntr)

        band_highlights = Dict(
            1:bandsbelow => attr(color=:red, width=3)
        )

        # return(plotbz(sgnum, Rs))

        p = plotbands(joinpath(out_dir, calcname), Rs, sg=sgnum, dim=3,
            band_highlights=band_highlights,
            title="$trialid-$localid-$id     Topological index: $topoindex      Compgap: $(@sprintf("%.4f", completegap))",
            savename=joinpath(out_dir_full, "bands-zoom-$id.svg"), nkpoints=300,
            range_fill=[gap_bounds[1], gap_bounds[2]], 
            ylim=[0.5, 0.63], 
            width=600, height=250
        )
        # display(p)

    end
end