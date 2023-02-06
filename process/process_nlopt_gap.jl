#=
Process NLopt results for optimizing bandgap.
=#
using MPBUtils
using DelimitedFiles
include("../src/bands.jl")
using FileIO
using Printf
include("../src/topology.jl")
include("../src/process_symeigs.jl")

sgnum = parse(Int64, ARGS[1])

n = 10000
# nin = 1.62
# nin = 3.4757
nin = 4

# sgnum = 13
D = 3

prefix = "opt-global-wideisoval-direct/trial"
saveprefix = "global"
mu1 = minconnectivity(sgnum)

# # SG 13
# trial_ids = [4, 5, 6, 7, 7, 7, 8, 8, 9]
# global_best_ids = [1106, 899, 856, 1670, 1973, 1673, 1397, 1677, 1950]
# # SG 68
# trial_ids = [1, 1, 1, 2, 2, 3, 3, 3, 4, 7, 8, 9]
# global_best_ids = [722, 772, 909, 1473, 1993, 800, 1111, 1016, 1610, 1150, 1351, 1836]

# if sgnum==68
#     # globalmin
#     # trial_ids = [1, 1, 1, 2, 2, 3, 3, 3, 4, 7, 8, 9]
#     # global_best_ids = [722, 772, 909, 1473, 1993, 800, 1111, 1016, 1610, 1150, 1351, 1836]
#     trial_ids = [2, 2, 3, 3, 3, 5, 6, 6, 7]
#     global_best_ids = [342, 1167, 1455, 1883, 1143, 1248, 1617, 844, 1724]
# elseif sgnum==13
#     trial_ids = [1, 2, 2, 3, 4, 4, 6, 7, 7, 9]
#     global_best_ids = [1560, 1582, 1736, 445, 579, 1587, 1804, 1242, 1572, 1916]
# elseif sgnum==48
#     trial_ids = [3, 4, 5, 5, 6, 6, 8]
#     global_best_ids = [1618, 1435, 1391, 1796, 956, 578, 160]
# elseif sgnum==49
#     trial_ids = [2, 2, 2, 3, 3, 8, 10]
#     global_best_ids = [805, 677, 107, 1860, 860, 933, 1151]
# elseif sgnum==50
#     trial_ids = [1, 2, 5, 5, 6, 7, 7, 8]
#     global_best_ids = [1777, 683, 1132, 1055, 1785, 1251, 1558, 534]
# elseif sgnum==86
#     trial_ids = [2, 2, 4, 4, 8, 8, 8, 8, 8, 9]
#     global_best_ids = [1541, 1466, 1645, 1528, 1307, 1545, 1845, 1457, 1169, 1532]
# end

for trial in 1:20
    # These are sub-folders in /input/ and /output/, respectively
    # Name is something like "sg68-n3.48-trial1"
    trialid = trial
    in_dir = "sg$sgnum-dim$D-n$(@sprintf("%.2f", nin))/$prefix$trialid"
    # trialid = trial_ids[trial]
    # localid = global_best_ids[trial]
    # in_dir = "sg$sgnum-dim$D-n$(@sprintf("%.2f", nin))/$prefix-$trialid-$localid"
    out_dir = in_dir
    in_dir_full = joinpath("input", in_dir)
    out_dir_full = joinpath("output", out_dir)
    if !isdir(in_dir_full)
        throw(ArgumentError("Directory $in_dir does not exist."))
    end
 
    #= Collect all the bandgaps
    Since each trial may have different number of iterations, we walk through files until they don't exist.
    File structure is different depending on whether the space group is filling-enforced or not.
    For non-filling-enforced space groups, there are many iterations without full dispersion calculations.
    =#
    objective_list = []
    objective2_list = []
    maxobj_list = []

    for i in 1:n
        calcname = "$i"
        if isfile(joinpath(out_dir_full, calcname*"-sym-dispersion.out"))
            if filesize(joinpath(out_dir_full, calcname*"-dispersion.out")) > 0
                gap_i = calccompletebandgap(joinpath(out_dir, calcname), mu1)
                gap2_i = calcminbandgap(joinpath(out_dir, calcname), mu1)
            else
                gap_i = -10
                gap2_i = -10
            end
            append!(objective_list, gap_i)
            append!(objective2_list, gap2_i)

            # Keep running track of best value so far
            if i==1 || gap_i > maxobj_list[end]
                append!(maxobj_list, gap_i)
            else
                append!(maxobj_list, maxobj_list[end])
            end
        else
            break
        end
    end

    println("Processed " * string(length(objective_list)) * " files.")
    println("Topological count: " * string(count(x -> x > -10, objective_list)))
    inds_max = partialsortperm(objective_list, 1:5, rev=true)     # Get top 5 bandgaps
    gaps_max = objective_list[inds_max]
    println(inds_max)
    println(gaps_max)
    println(objective2_list[inds_max])

    # Write results of top 5 candidates to a file
    out_dir_top = splitdir(out_dir_full)[1]
    open(joinpath(out_dir_top, "$saveprefix-summary.txt"), "a") do io
        println(io, "Trial $trialid, in_dir $in_dir")
        println(io, "Processed " * string(length(objective_list)) * " files.")
        println(io, "Topological count: " * string(count(x -> x > -10, objective_list)))
        println(io, inds_max)
        println(io, gaps_max)
        println(io, gaps_max)
        println(io, objective2_list[inds_max])
    end

    open(joinpath(out_dir_top, "$saveprefix-obj-best.txt"), "a") do io
        writedlm(io, transpose(maxobj_list))
    end

    for id in inds_max
        if objective_list[id] > -10
            calcname = "$id"
            best_dir = in_dir*"-best"
            in_dir_best_full = joinpath("input", best_dir)
            out_dir_best_full = joinpath("output", best_dir)
            isdir(in_dir_best_full ) || mkpath(in_dir_best_full)
            isdir(out_dir_best_full) || mkpath(out_dir_best_full)
            cp(joinpath(in_dir_full, calcname*"-sym.sh"),
                joinpath(in_dir_best_full, calcname*"-sym.sh"), force=true)
            cp(joinpath(in_dir_full, calcname*".sh"),
                joinpath(in_dir_best_full, calcname*".sh"), force=true)
            cp(joinpath(out_dir_full, calcname*"-sym-symeigs.out"),
                joinpath(out_dir_best_full, calcname*"-sym-symeigs.out"), force=true)
            cp(joinpath(out_dir_full, calcname*"-sym-dispersion.out"),
                joinpath(out_dir_best_full, calcname*"-sym-dispersion.out"), force=true)
            cp(joinpath(out_dir_full, calcname*"-dispersion.out"),
                joinpath(out_dir_best_full, calcname*"-dispersion.out"), force=true)
        end
    end
end