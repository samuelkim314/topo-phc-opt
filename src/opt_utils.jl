# Collection of methods to facilitate optimization, track candidates, and interface with NLopt.
using DelimitedFiles
using Printf
using NLopt: LN_SBPLX
using NLopt: GN_ISRES
using NLopt: GN_DIRECT_L_RAND
using NLopt

mutable struct OptTracker
    x_best::Vector
    id_best::Integer
    obj_best::Real
    x_list::Vector
    id_list::Vector
    obj_list::Vector
    x_best_list::Vector
    id_best_list::Vector
    obj_best_list::Vector
    x_topk::Vector
    id_topk::Vector
    obj_topk::Vector
end
OptTracker() = OptTracker([], 0, -2.0, [], [], [], [], [], [], [], [], [])


"""
Update the `tracker` object with the most recent evaluation so it can update the list of best candidate(s).

# Arguments
- `x::Vector`: Parameter vector of the candidate
- `id::Integer`: ID number corresponding to the candidate
- `obj::Real`: Objective value (e.g. bandgap)
- `constraint::Bool`: whether the candidate obeys the constraint
- `nkeep::Integer`: number of candidates to keep track of
"""
function adddata!(tracker::OptTracker, x::Vector, id::Integer, obj::Real, constraint::Bool=true;
    nkeep::Integer=5)
    x = copy(x) # x is mutable list, so will just make a pointer otherwise
    
    push!(tracker.x_list, x)
    append!(tracker.id_list, id)
    append!(tracker.obj_list, obj)

    if constraint
        if obj > tracker.obj_best
            tracker.x_best = x
            tracker.obj_best = obj
            tracker.id_best = id
        end

        if length(tracker.obj_topk) < nkeep
            if obj > -2.0
                push!(tracker.x_topk, x)
                append!(tracker.obj_topk, obj)
                append!(tracker.id_topk, id)
            end
        else
            jmin = argmin(tracker.obj_topk)
            if obj > tracker.obj_topk[jmin]
                tracker.x_topk[jmin] = x
                tracker.obj_topk[jmin] = obj
                tracker.id_topk[jmin] = id
            end
        end
        
    end
    # Keep a running list of the best value found so far
    push!(tracker.x_best_list, tracker.x_best)
    append!(tracker.obj_best_list, tracker.obj_best)
    append!(tracker.id_best_list, tracker.id_best)
end


"""
    recordbest(tracker, indir, outdir, outdirbest[, rmold])

Move the input/output files corresponding to the top candidates stored in `tracker` from `indir`/`outdir` to 
`outdirbest`.

If `rmold` is `true`, then it will clear the `outdirbest` directory first before moving the candidate files.
This is convenient to avoid piling up hundreds of unnecessary files.
"""
function recordbest(tracker::OptTracker, indir::String, outdir::String, outdirbest::String=nothing;
    rmold::Bool=true)

    if isnothing(outdirbest)
        outdirbest = outdir*"-best"
    end

    if rmold
        # Remove all the existing .sh and .out files to keep it clean
        fileslist = readdir(outdirbest, join=true)
        for file in fileslist
            if endswith(file, ".out") || endswith(file, ".sh")
                rm(file)
            end
        end
    end

    # Record best values so far
    open(joinpath(outdirbest, "summary.txt"), "w") do io
        println(io, "out_dir $outdirbest")
        println(io, "Ran " * string(length(tracker.obj_list)) * " iterations.")
        println(io, "Topological count: " * string(count(x -> x > -2, tracker.obj_list)))
        println(io, tracker.id_topk)
        println(io, tracker.obj_topk)
    end
    # Record running list of best value found so far - makes it convenient for plotting
    open(joinpath(outdirbest, "obj-best.txt"), "w") do io
        writedlm(io, transpose(tracker.obj_best_list))
    end
    # Record full data
    open(joinpath(outdirbest, "fulldatax.txt"), "w") do io
        writedlm(io, tracker.x_list)
    end
    open(joinpath(outdirbest, "fulldataobj.txt"), "w") do io
        writedlm(io, tracker.obj_list)
    end

    # Copy over all the MPB simulation and output files
    for j in tracker.id_topk
        cp(joinpath(indir, "$j-sym.sh"),
            joinpath(outdirbest, "$j-sym.sh"), force=true)
        cp(joinpath(indir, "$j.sh"),
            joinpath(outdirbest, "$j.sh"), force=true)
        cp(joinpath(outdir, "$j-sym-symeigs.out"),
            joinpath(outdirbest, "$j-sym-symeigs.out"), force=true)
        cp(joinpath(outdir, "$j-sym-dispersion.out"),
            joinpath(outdirbest, "$j-sym-dispersion.out"), force=true)
        cp(joinpath(outdir, "$j-dispersion.out"),
            joinpath(outdirbest, "$j-dispersion.out"), force=true)
    end
end


"""
    in_dir, out_dir, log_dir, in_dir_full, out_dir_full, out_dir_best_full = makedirs(subdir[, tmpdir])

Create directories for saving input/output/log files. 

`in_dir` contains .sh files for specifying PhC parameters.
`out_dir` contains output files, including dispersion files and symmetry eigenvalues.
`log_dir` contains MPB output (from which dispersion and symmetry eigenvalues are extracted)
If all of `in_dir`, `out_dir`, and `log_dir` are provided (as keyword arguments), then it will use these.
Otherwise ff `tmpdir == true` and environment variable TMPDIR is available, it will use TMPDIR for input/output/log
directories.
This is convenient for cluster configurations where TMPDIR is a local directory not synchronized with the cluster.
Otherwise, it uses "input/", "output/", and "logs/" by default.

`subdir` is typically something like "sg68-dim3-n4.00/chern/trial1"

`out_dir_best_full` does not use TMPDIR and uses "output" since we want to save these results.
"""
function makedirs(subdir::String, tmpdir::Bool=false; 
    in_dir=nothing, out_dir=nothing, log_dir=nothing, bestsuffix::String="-best")

    if isnothing(in_dir) || isnothing(out_dir) || isnothing(log_dir)
        if tmpdir && haskey(ENV, "TMPDIR")
            in_dir = ENV["TMPDIR"]
            out_dir = ENV["TMPDIR"]
            log_dir = ENV["TMPDIR"]
        else
            in_dir = "input"
            out_dir = "output"
            log_dir = "logs"
        end
    end
    
    in_dir_full = joinpath(in_dir, subdir)
    out_dir_full = joinpath(out_dir, subdir)
    out_dir_best_full = joinpath("output", subdir*bestsuffix)
    isdir(in_dir_full ) || mkpath(in_dir_full)
    isdir(out_dir_full) || mkpath(out_dir_full)
    isdir(out_dir_best_full) || mkpath(out_dir_best_full)
    isdir(joinpath(log_dir, subdir)) || mkpath(joinpath(log_dir, subdir))

    return in_dir, out_dir, log_dir, in_dir_full, out_dir_full, out_dir_best_full
end


"""
Global optimization with proper bounds for photonic crystal lattice-set parameterization.

# Arguments
- `f::Function`: function to optimize with signature `f(x::Vector, grad::Vector, opt_id::Integer)` and returns the objective.
    `globalopt` puts `f` inside a wrapper to track the best candidates and throw any errors 
    (since NLopt suppresses errors).
- `x0::Vector`: Initial candidate
- `nlatticeparams::Integer`: Number of lattice parameters, see `nbasisparams` in `params.jl`
- `nflatparams::Integer`: Number of free Fourier sum coefficients
- `hasinv::Bool`: Has inversion symmetry.
- `maxeval::Integer`: Number of evaluations for global optimization.
"""
function globalopt_levelset(f::Function, x0::Vector, nlatticeparams::Integer, nflatparams::Integer, hasinv::Bool, 
    maxeval::Integer, in_dir_full, out_dir_full, out_dir_best_full)
    
    opt_id = 1
    tracker = OptTracker()

    """
    Wrapper for the function be optimized that trackers the best candidates.
    """
    function ftracked(x::Vector, grad::Vector)
        println("Iteration $opt_id")
        try
            bandgap = f(x, grad, opt_id)

            # Collect running statistics on best data
            adddata!(tracker, x, opt_id, bandgap)
            println("Iteration $opt_id finished.\nObjective : $bandgap\nx: $x")

            if opt_id % 500 == 0
                recordbest(tracker, in_dir_full, out_dir_full, out_dir_best_full)
            end

            opt_id += 1
            return bandgap
        catch e
            # Throw the error manually because NLopt will suppress the error otherwise
            println(e)
            @error "Something went wrong" exception=(e, catch_backtrace())
            throw(e)
        end
    end

    opt = Opt(GN_DIRECT_L_RAND, length(x0))
    # opt.population = 50
    opt.maxeval = maxeval
    if hasinv
        opt.lower_bounds = append!([0], zeros(nlatticeparams), -1*ones(nflatparams))
        opt.upper_bounds = append!([1], ones(nlatticeparams), ones(nflatparams))
    else
        opt.lower_bounds = append!([0], zeros(nlatticeparams), -1*ones(2*nflatparams))
        opt.upper_bounds = append!([1], ones(nlatticeparams), ones(2*nflatparams))
    end

    opt.max_objective = ftracked

    (optf,optx,ret) = optimize(opt, x0)

    println(optf)
    println(optx)
    println(ret)

    return tracker
end


"""
Stochastic local optimization using a combination of SBPLX and ISRES.

# Arguments
- `f`: function to optimize, in accordance to NLopt format. Takes x and gradient as input, and returns objective.
- `nlatticeparams::Integer`: Number of lattice parameters
- `nflatparams::Integer`: Number of free Fourier coefficients
- `hasinv::Bool`: whether the structure has inversion symmetry
- `ncycles::Integer`: number of optimization cycles
- `iters_per_cycle::Integer`: number of optimization iterations per cycle
"""
function localopt_levelset(f, x0::Vector, nlatticeparams::Integer, nflatparams::Integer, hasinv::Bool,
    ncycles::Integer, iters_per_cycle::Integer,
    in_dir_full::String, out_dir_full::String, out_dir_best_full::String)

    opt_id = 1
    tracker = OptTracker()

    """
    Wrapper for the function be optimized that trackers the best candidates.
    """
    function ftracked(x::Vector, grad::Vector)
        println("Iteration $opt_id")
        try
            bandgap = f(x, grad, opt_id)

            # Collect running statistics on best data
            adddata!(tracker, x, opt_id, bandgap)
            println("Iteration $opt_id finished.\nObjective : $bandgap\nx: $x")

            opt_id += 1
            return bandgap
        catch e
            # Throw the error manually because NLopt will suppress the error otherwise
            println(e)
            @error "Something went wrong" exception=(e, catch_backtrace())
            throw(e)
        end
    end

    # Set up optimizer parameters
    opt1 = Opt(LN_SBPLX, length(x0))
    opt2 = Opt(GN_ISRES, length(x0))

    opt1.maxeval = iters_per_cycle
    opt2.maxeval = iters_per_cycle

    # Optimization bounds. 
    if hasinv
        # Inversion symmetry means that coefficients are restricted to be real
        lower_bounds = append!([0], zeros(nlatticeparams), -1*ones(nflatparams))
        upper_bounds = append!([1], ones(nlatticeparams), ones(nflatparams))
    else
        # Fourier coefficients are, so we split into real + imaginary components.
        lower_bounds = append!([0], zeros(nlatticeparams), -1*ones(2*nflatparams))
        upper_bounds = append!([1], ones(nlatticeparams), ones(2*nflatparams))
    end

    # Setting initial step size
    range = upper_bounds - lower_bounds
    range2 = copy(range)
    range2[1] = 0   # We don't want to change bounds for ff
    opt1.initial_step = range / 10

    opt1.max_objective = ftracked
    opt2.max_objective = ftracked

    for i in 1:ncycles
        id_best = tracker.id_best
        obj_best = tracker.obj_best
        println("Best before cycle $i: $id_best, $obj_best")
        println("Starting point: $x0")

        if id_best <= (i-3) * iters_per_cycle # if previous 2 cycles saw no improvement
            println("Entering wide global search")
            opt2.lower_bounds = x0 - range / 10
            opt2.upper_bounds = x0 + range / 10
            opt2.lower_bounds[1] = 0
            opt2.upper_bounds[1] = 1
            (optf,optx,ret) = optimize(opt2, x0)
            upper_bounds += range2/10
            lower_bounds -= range2/10
        elseif id_best <= (i-2) * iters_per_cycle # if previous cycle saw no improvement
            println("Entering global search")
            opt2.lower_bounds = x0 - range / 20
            opt2.upper_bounds = x0 + range / 20
            opt2.lower_bounds[1] = 0
            opt2.upper_bounds[1] = 1
            (optf,optx,ret) = optimize(opt2, x0)
            upper_bounds += range2/20
            lower_bounds -= range2/20
        else
            println("Local search")
            # Set bounds for the local optimizer
            opt1.lower_bounds = lower_bounds
            opt1.upper_bounds = upper_bounds
            (optf,optx,ret) = optimize(opt1, x0)
        end

        id_best = tracker.id_best
        obj_best = tracker.obj_best
        println("Best after cycle $i: $id_best, $obj_best")
        println(tracker.x_best)

        println(optf)
        println(optx)
        println(ret)

        # Set starting point for next cycle
        if !isempty(tracker.x_best)
            x0 = copy(tracker.x_best)
        end

        recordbest(tracker, in_dir_full, out_dir_full, out_dir_best_full)
    end
end


"""
Global optimization with [0, 1] bounds.

# Arguments
- `f::Function`: function to optimize with signature `f(x::Vector, grad::Vector, opt_id::Integer)` and returns the objective.
    `globalopt` puts `f` inside a wrapper to track the best candidates and throw any errors 
    (since NLopt suppresses errors).
- `x0::Vector`: Initial candidate
- `nlatticeparams::Integer`: Number of lattice parameters, see `nbasisparams` in `params.jl`
- `nflatparams::Integer`: Number of free Fourier sum coefficients
- `hasinv::Bool`: Has inversion symmetry.
- `maxeval::Integer`: Number of evaluations for global optimization.
"""
function globalopt(f::Function, x0::Vector, 
    maxeval::Integer, in_dir_full, out_dir_full, out_dir_best_full)
    
    opt_id = 1
    tracker = OptTracker()

    """
    Wrapper for the function be optimized that trackers the best candidates.
    """
    function ftracked(x::Vector, grad::Vector)
        println("Iteration $opt_id")
        try
            bandgap = f(x, grad, opt_id)

            # Collect running statistics on best data
            adddata!(tracker, x, opt_id, bandgap)
            println("Iteration $opt_id finished.\nObjective : $bandgap\nx: $x")

            if opt_id % 500 == 0
                recordbest(tracker, in_dir_full, out_dir_full, out_dir_best_full)
            end

            opt_id += 1
            return bandgap
        catch e
            # Throw the error manually because NLopt will suppress the error otherwise
            println(e)
            @error "Something went wrong" exception=(e, catch_backtrace())
            throw(e)
        end
    end

    opt = Opt(GN_DIRECT_L_RAND, length(x0))
    # opt.population = 50
    opt.maxeval = maxeval
    nparams = length(x0)
    opt.lower_bounds = zeros(nparams)
    opt.upper_bounds = ones(nparams)

    opt.max_objective = ftracked

    (optf,optx,ret) = optimize(opt, x0)

    println(optf)
    println(optx)
    println(ret)

    return tracker
end


"""
Stochastic local optimization using a combination of SBPLX and ISRES.

# Arguments
- `f`: function to optimize, in accordance to NLopt format. Takes x and gradient as input, and returns objective.
- `nlatticeparams::Integer`: Number of lattice parameters
- `nflatparams::Integer`: Number of free Fourier coefficients
- `hasinv::Bool`: whether the structure has inversion symmetry
- `ncycles::Integer`: number of optimization cycles
- `iters_per_cycle::Integer`: number of optimization iterations per cycle
"""
function localopt(f, x0::Vector, 
    ncycles::Integer, iters_per_cycle::Integer,
    in_dir_full::String, out_dir_full::String, out_dir_best_full::String)

    opt_id = 1
    tracker = OptTracker()

    """
    Wrapper for the function be optimized that trackers the best candidates.
    """
    function ftracked(x::Vector, grad::Vector)
        println("Iteration $opt_id")
        try
            bandgap = f(x, grad, opt_id)

            # Collect running statistics on best data
            adddata!(tracker, x, opt_id, bandgap)
            println("Iteration $opt_id finished.\nObjective : $bandgap\nx: $x")

            opt_id += 1
            return bandgap
        catch e
            # Throw the error manually because NLopt will suppress the error otherwise
            println(e)
            @error "Something went wrong" exception=(e, catch_backtrace())
            throw(e)
        end
    end

    # Set up optimizer parameters
    opt1 = Opt(LN_SBPLX, length(x0))
    opt2 = Opt(GN_ISRES, length(x0))

    opt1.maxeval = iters_per_cycle
    opt2.maxeval = iters_per_cycle

    # "Global" Optimization bounds. 
    nparams = length(x0)
    lower_bounds = zeros(nparams)
    upper_bounds = ones(nparams)
    # Set bounds for the local optimizer
    opt1.lower_bounds = lower_bounds
    opt1.upper_bounds = upper_bounds

    # Setting initial step size
    range = upper_bounds - lower_bounds
    opt1.initial_step = range / 10

    opt1.max_objective = ftracked
    opt2.max_objective = ftracked

    for i in 1:ncycles
        id_best = tracker.id_best
        obj_best = tracker.obj_best
        println("Best before cycle $i: $id_best, $obj_best")
        println("Starting point: $x0")

        if id_best <= (i-3) * iters_per_cycle # if previous 2 cycles saw no improvement
            println("Entering wide global search")
            opt2.lower_bounds = max.(upper_bounds, x0 - range / 10)
            opt2.upper_bounds = min.(lower_bounds, x0 + range / 10)
            (optf,optx,ret) = optimize(opt2, x0)
        elseif id_best <= (i-2) * iters_per_cycle # if previous cycle saw no improvement
            println("Entering global search")
            opt2.lower_bounds = max.(upper_bounds, x0 - range / 20)
            opt2.upper_bounds = min.(lower_bounds, x0 + range / 20)
            (optf,optx,ret) = optimize(opt2, x0)
        else
            println("Local search")
            (optf,optx,ret) = optimize(opt1, x0)
        end

        id_best = tracker.id_best
        obj_best = tracker.obj_best
        println("Best after cycle $i: $id_best, $obj_best")
        println(tracker.x_best)

        println(optf)
        println(optx)
        println(ret)

        # Set starting point for next cycle
        if !isempty(tracker.x_best)
            x0 = copy(tracker.x_best)
        end

        recordbest(tracker, in_dir_full, out_dir_full, out_dir_best_full)
    end
end