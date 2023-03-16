#=
Methods for getting bands from MPB calculations, calculating objective metrics, and plotting
    band structures.
=#

using Crystalline
using MPBUtils
using Brillouin
using DelimitedFiles
using PlotlyJS
using FileIO
include("topology.jl")
include("write_ctl.jl")


function getbands(calcname::String; dir::String="output")
    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"),',')
    freqs = dispersion_data[:,6:end]
    return freqs
end

function getband(calcname::String, bandnum::Union{Integer, AbstractArray};
        dir::String="output")
    return getbands(calcname; dir=dir)[:, bandnum]
end

function getgap(calcname::String, lower::Integer; dir::String="output")
    freqs = getbands(calcname; dir=dir)
    gap_min = maximum(freqs[:, 1:lower])    # bottom of the gap
    gap_max = minimum(freqs[:, lower+1:end])    # top of the gap
    return [gap_min, gap_max]
end


"""
    calcminbandgap(calcname[, lower]; dir)

Calculate minimum pointwise bandgap for MPB output file `calcname`.

If `lower` is an Integer, compute the min bandgap above band `lower`.
If `lower` is a list of UnitRanges (that typically represent connected bands), then
compute the min bandgap above all the UnitRanges and return the largest one along with the 
index of the UnitRange in the `lower` list corresponding to that gap.
If `lower` is unspecified, compute bandgap between all pairs of adjacent bands and 
return the largest one along with the index of the band below the gap.
"""
function calcminbandgap(freqs::Matrix, lower::Integer=2)
    if lower <= 0 || lower >= size(freqs, 2)
        throw(ArgumentError("lower must be between 1 and nbands-1, 
            where nbands=size(freqs, 2)"))
    end

    gap_min = maximum(freqs[:, 1:lower], dims=2)    # bottom of the gap array
    gap_max = minimum(freqs[:, lower+1:end], dims=2)    # top of the gap array

    # Bandgap calculated as a fraction of mid-gap frequency
    bandgap = (gap_max .- gap_min) ./ (gap_max .+ gap_min) .* 2
    minbandgap = minimum(vec(bandgap))

    return minbandgap
end
function calcminbandgap(calcname::String, lower::Integer=2; dir::String="output")
    freqs = getbands(calcname; dir=dir)
    return calcminbandgap(freqs, lower)
end
function calcminbandgap(calcname::String, lowerbands::Vector; dir::String="output")
    freqs = getbands(calcname; dir=dir)

    obj_list = []
    for bands in lowerbands
        if bands[end] >= size(freqs, 2)
            append!(obj_list, -2.0)
        else
            append!(obj_list, calcminbandgap(freqs, bands[end]))
        end
    end
    return findmax(obj_list)
end
function calcminbandgap(calcname::String; dir::String="output")
    freqs = getbands(calcname; dir=dir)
    n_bands = size(freqs, 1)
    obj_list = []
    for lower in 1:n_bands-1
        append!(obj_list, calcminbandgap(freqs, lower))
    end
    return findmax(obj_list)
end

"""
Find the location along k with the min bandgap
"""
function calcargminbandgap(calcname, lower::Integer=2)
    dispersion_data = readdlm((@__DIR__)*"/output/"*calcname*"-dispersion.out", ',')
    freqs = dispersion_data[:,6:end]

    gap_min = nothing
    gap_min = maximum(freqs[:, 1:lower], dims=2)    # bottom of the gap array
    gap_max = minimum(freqs[:, lower+1:end], dims=2)    # top of the gap array

    # Bandgap calculated as a fraction of mid-gap frequency
    bandgap = (gap_max .- gap_min) ./ (gap_max .+ gap_min) .* 2
    minbandgap, i = findmin(vec(bandgap))

    return minbandgap, i
end

"""
    calccompletebandgap(freqs[, lower])
    calccompletebandgap(calcname[, lower])

Calculate complete bandgap for MPB output file `calcname` or matrix `freqs`.
Looks at the band between lower and lower+1. lower must be between 1 and size(freqs,2)-1
If `lower` is a vector of UnitRanges, than take the maximum gap of all the ends of ranges in `lower`
"""
function calccompletebandgap(freqs::Matrix, lower::Integer=2)
    if lower <= 0 || lower >= size(freqs, 2)
        throw(ArgumentError("lower must be between 1 and nbands-1, 
            where nbands=size(freqs, 2)"))
    end

    max_lower = maximum(freqs[:, 1:lower])  # Max of nth band
    min_upper = minimum(freqs[:, lower+1:end])  # Min of (n+1)th band

    # Bandgap calculated as a fraction of mid-gap frequency
    bandgap = (min_upper - max_lower) / (min_upper + max_lower) * 2

    return bandgap
end
function calccompletebandgap(calcname::String, lower::Integer=2; dir::String="output")
    freqs = getbands(calcname; dir=dir)
    return calccompletebandgap(freqs, lower)
end
function calccompletebandgap(calcname::String, lowerbands::Vector; dir::String="output")
    freqs = getbands(calcname; dir=dir)

    obj_list = []
    for bands in lowerbands
        if bands[end] >= size(freqs, 2)
            append!(obj_list, -2.0)
        else
            append!(obj_list, calccompletebandgap(freqs, bands[end]))
        end
    end
    return findmax(obj_list)
end


"""
Min gap metric but ignoring k-paths adjacent to the `label_skip` high-symmetry k-point.
Useful for optimization when a degeneracy lies on the high-symmetry k-path.
"""
function partialmingap(kpi, freqs, label_skip; degen_bands=1:6, nbands_below=3)
    local_xs         = cumdists.(kpi.kpaths)  
    start_idx = 1
    degen_idx = Int32[]
    nondegen_idx = Int32[]
    for (path_idx, (local_x, labels)) in enumerate(zip(local_xs, kpi.labels))
        stop_idx = start_idx+length(local_x)-1
        append!(degen_idx, [k+start_idx-1 for (k, v) in labels if v==label_skip])
        append!(nondegen_idx, [k+start_idx-1 for (k, v) in labels if v!=label_skip])
        start_idx += stop_idx
    end

    n_bands = size(freqs, 2)

    gap_bottom_idx = degen_bands[1] + nbands_below - 1
    gap_top_idx = gap_bottom_idx + 1

    gap_bottom = maximum(freqs[nondegen_idx, 1:gap_bottom_idx], dims=2)
    gap_top = minimum(freqs[nondegen_idx, gap_top_idx:size(freqs, 2)], dims=2)

    w0 = freqs[degen_idx[1], degen_bands[1]]

    gap = gap_top - gap_bottom
    mingap = minimum(gap) / w0

    return mingap
end

"""
Complete gap metric but ignoring k-paths adjacent to the `label_skip` high-symmetry k-point.
Useful for optimization when a degeneracy lies on the high-symmetry k-path.
"""
function partialcompletegap(kpi, freqs, label_skip; degen_bands=1:6, nbands_below=3)
    local_xs         = cumdists.(kpi.kpaths)  
    start_idx = 1
    degen_idx = Int32[]
    nondegen_idx = Int32[]
    for (path_idx, (local_x, labels)) in enumerate(zip(local_xs, kpi.labels))
        stop_idx = start_idx+length(local_x)-1
        append!(degen_idx, [k+start_idx-1 for (k, v) in labels if v==label_skip])
        append!(nondegen_idx, [k+start_idx-1 for (k, v) in labels if v!=label_skip])
        start_idx += stop_idx
    end

    n_bands = size(freqs, 2)

    gap_bottom_idx = degen_bands[1] + nbands_below - 1
    gap_top_idx = gap_bottom_idx + 1

    gap_bottom = maximum(freqs[nondegen_idx, 1:gap_bottom_idx])
    gap_top = minimum(freqs[nondegen_idx, gap_top_idx:n_bands])
    
    w0 = freqs[degen_idx[1], degen_bands[1]]

    completegap = (gap_top - gap_bottom) / w0

    return completegap
end
function partialcompletegapanyband(kpi, freqs, label_skip; degen_bands=1:6)
    gap_list = [partialcompletegap(kpi, freqs, label_skip, degen_bands=degen_bands, nbands_below=i) 
        for i in 1:5]
    gapmax, imax = findmax(gap_list)

    return gapmax, imax
end


"""
    calcpartialhinge(kpi, freqs, degenbands, klabskip)

Calculate partial hinge bandgap metric
"""
function calcpartialhinge(kpi::KPathInterpolant{3}, freqs::Matrix, 
    degenbands::Union{Vector{Integer}, UnitRange}, klabskip::Symbol,
    nbandsbelow::Integer=3)

    # Get array indices corresponding to the k-points we want to calculate at
    local_xs = cumdists.(kpi.kpaths)  
    start_idx = 1
    degen_idx = Int32[]
    nondegen_idx = Int32[]
    for (local_x, labels) in zip(local_xs, kpi.labels)
        stop_idx = start_idx+length(local_x)-1
        append!(degen_idx, [k+start_idx-1 for (k, v) in labels if v==klabskip])
        append!(nondegen_idx, [k+start_idx-1 for (k, v) in labels if v!=klabskip])
        start_idx = stop_idx + 1
    end

    nondegen_idx = deleteat!(collect(1:length(kpi)), degen_idx)

    # Find the top and bottom of the "gap"
    n_bands = size(freqs, 2)
    gap_bottom_idx = degenbands[1] + nbandsbelow - 1
    gap_top_idx = gap_bottom_idx + 1
    gap_bottom = maximum(freqs[nondegen_idx, 1:gap_bottom_idx])
    gap_top = minimum(freqs[nondegen_idx, gap_top_idx:n_bands])

    w0 = freqs[degen_idx[1], degenbands[1]]    # Frequency of the degeneracy

    # Applying the objective function
    abovehinge = min(0, (gap_top - w0) / w0)
    belowhinge = min(0, (w0 - gap_bottom) / w0)
    gap = abovehinge + belowhinge

    return gap
end

function calcpartialhinge(kpi::KPathInterpolant{3}, freqs::Matrix, 
    degenbands::Union{Vector{Integer}, UnitRange}, klabskip::Symbol)
    degendim = maximum(degenbands) - minimum(degenbands) + 1

    gap_list = [calcpartialhinge(kpi, freqs, degenbands, klabskip, i) 
        for i in 1:degendim-1]
    gapmax, imax = findmax(gap_list)

    return gapmax, imax
end

function calcpartialhinge(kpi::KPathInterpolant{3}, freqs::Matrix, 
    degenbands::Vector{UnitRange}, klabskip::Symbol)

    gapband_list = [calcpartialhinge(kpi, freqs, degenbandsi, klabskip) 
        for degenbandsi in degenbands]
    # In Julia 1.7 this can be simplified, but findmax has different behavior in 1.6
    gap_list = [i[1] for i in gapband_list]
    gapmax, imax = findmax(gap_list)
    
    return gapmax, gapband_list[imax][2], degenbands[imax]
end

"""
TODO: Finish testing this function (mostly the k-point interpolation)
"""
function calcpartialhinge(incalcpath::String, freqs::Matrix, degenbands::Vector, 
    sgnum::Integer=nothing, klabskip::Symbol=nothing; nkpoints=nothing)

    if isnothing(sgnum)
        sgnum = try_parse_sgnum(calcname)    # int for space group
    end
    if isnothing(klabskip)
        klabskip = getdegenksym(sgnum)
    end
    if isnothing(nkpoints)
        nkpoints = size(freqs, 2)
    end

    pRs, _, _, _, _, _ = lattice_from_mpbparams(incalcpath)
    cntr = centering(sgnum, 3) # will be ‘C’ for base-centered
    Rs = conventionalize(pRs, cntr)
    kp = irrfbz_path(sgnum, Rs)
    kpi = interpolate(kp, nkpoints)
    # kpi = kp

    return calcpartialhinge(kpi, freqs, degenbands, klabskip)

end

"""
    calcgaphinge(freqs, nbandsbelow, gap0lower, gap0upper)
    calcgaphinge(calcname, nbandsbelow, gap0lower, gap0upper)

Calculate hinge loss for a complete bandgap in a specified region. Useful for
optimizing a cladding material.
"""
function calcgaphinge(freqs::Matrix, nbandsbelow::Integer, gap0lower::Real, gap0upper::Real)

    gap_lower = maximum(freqs[:, 1:nbandsbelow], dims=2)
    gap_upper = minimum(freqs[:, nbandsbelow+1:end], dims=2)
    max_lower = maximum(gap_lower)
    min_upper = minimum(gap_upper)

    # I don't think we need to normalize by w0 since gap0 is fixed
    abovehinge = min(0, (min_upper - gap0upper))
    belowhinge = min(0, (gap0lower - max_lower))
    objective = abovehinge + belowhinge

    return objective
end
function calcgaphinge(calcname::String, nbandsbelow::Integer, gap0lower::Real, gap0upper::Real;
        dir::String="output")
    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"),',')
    freqs = dispersion_data[:,6:end]

    return calcgaphinge(freqs, nbandsbelow, gap0lower, gap0upper)
end
function calcgaphinge(calcname::String, nbandsbelow::Vector, gap0lower::Real, gap0upper::Real;
        dir::String="output")
    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"),',')
    freqs = dispersion_data[:,6:end]

    obj_list = []
    for bands in nbandsbelow
        append!(obj_list, calcgaphinge(freqs, bands[end], gap0lower, gap0upper))
    end
    return findmax(obj_list)
end

"""
Calculate constraint that a bandgap completely encapsulates bounds `gap0lower` and `gap0upper`
For constraints of the form f_c(x)<=0
Useful for optimizing a cladding material
"""
function calcgapconstraint(freqs::Matrix, nbandsbelow::Integer, gap0lower::Real, gap0upper::Real)
    gap_lower = maximum(freqs[:, 1:nbandsbelow], dims=2)
    gap_upper = minimum(freqs[:, nbandsbelow+1:end], dims=2)
    max_lower = maximum(gap_lower)
    min_upper = minimum(gap_upper)

    return [max_lower - gap0lower, gap0upper - min_upper]
end
function calcgapconstraint(calcname::String, nbandsbelow::Integer, gap0lower::Real, gap0upper::Real;
        dir="output")
    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"),',')
    freqs = dispersion_data[:,6:end]

    return calcgapconstraint(freqs, nbandsbelow, gap0lower, gap0upper)
end

"""
Helper function that returns true if all constrains are fulfilled
"""
function allconstraints(calcname::String, nbandsbelow::Integer, gap0lower::Real, gap0upper::Real;
    dir="output")
    return all(calcgapconstraint(calcname, nbandsbelow, gap0lower, gap0upper; dir=dir) .<= 0)
end

"""
Calculate complete gap only if it fulfills constraints. Otherwise return -2.
"""
function calcconstrainedcompletegap(freqs::Matrix, nbandsbelow::Integer, gap0lower::Real, gap0upper::Real)
    constraints = calcgapconstraint(freqs, nbandsbelow, gap0lower, gap0upper)
    if all(constraints .<= 0)
        return calccompletebandgap(freqs, nbandsbelow)
    else
        return -2.0
    end
end
function calcconstrainedcompletegap(calcname::String, nbandsbelow::Integer, gap0lower::Real, gap0upper::Real;
        dir="output")
    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"),',')
    freqs = dispersion_data[:,6:end]

    return calcconstrainedcompletegap(freqs, nbandsbelow, gap0lower, gap0upper)
end
function calcconstrainedcompletegap(calcname::String, nbandsbelow::Vector, gap0lower::Real, gap0upper::Real;
    dir="output")
    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"),',')
    freqs = dispersion_data[:,6:end]

    obj_list = []
    for bands in nbandsbelow
        append!(obj_list, calcconstrainedcompletegap(freqs, bands[end], gap0lower, gap0upper))
    end
    return findmax(obj_list)
end

"""
Calculate complete gap only if it fulfills constraints. Otherwise return hinge loss
"""
function calchybridhingegap(freqs::Matrix, nbandsbelow::Integer, gap0lower::Real, gap0upper::Real)
    constraints = calcgapconstraint(freqs, nbandsbelow, gap0lower, gap0upper)
    if all(constraints .<= 0)
        return calccompletebandgap(freqs, nbandsbelow)
    else
        return calcgaphinge(freqs, nbandsbelow, gap0lower, gap0upper)
    end
end
function calchybridhingegap(calcname::String, nbandsbelow::Integer, gap0lower::Real, gap0upper::Real;
        dir="output")
    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"),',')
    freqs = dispersion_data[:,6:end]

    return calchybridhingegap(freqs, nbandsbelow, gap0lower, gap0upper)
end
function calchybridhingegap(calcname::String, nbandsbelow::Vector, gap0lower::Real, gap0upper::Real;
    dir="output")
    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"),',')
    freqs = dispersion_data[:,6:end]

    obj_list = []
    for bands in nbandsbelow
        append!(obj_list, calchybridhingegap(freqs, bands[end], gap0lower, gap0upper))
    end
    return findmax(obj_list)
end

"""
Calculate min cladding
"""
function calcmincladding(freqs::Matrix, nbandsbelow::Integer, gap0lower::Real, gap0upper::Real)
    if nbandsbelow <= 0 || nbandsbelow >= size(freqs, 2)
        throw(ArgumentError("nbandsbelow must be between 1 and nbands-1, 
            where nbands=size(freqs, 2)"))
    end

    gap_lower = maximum(freqs[:, 1:nbandsbelow], dims=2)
    gap_upper = minimum(freqs[:, nbandsbelow+1:end], dims=2)
    max_lower = maximum(gap_lower)
    min_upper = minimum(gap_upper)

    abovecladding = min_upper - gap0upper
    belowcladding = gap0lower - max_lower

    return min(abovecladding, belowcladding)
end
function calcmincladding(calcname::String, nbandsbelow::Integer, gap0lower::Real, gap0upper::Real;
    dir="output")

    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"),',')
    freqs = dispersion_data[:,6:end]

    return calcmincladding(freqs, nbandsbelow, gap0lower, gap0upper)
end
function calcmincladding(calcname::String, nbandsbelow::Vector, gap0lower::Real, gap0upper::Real;
    dir="output")
    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"),',')
    freqs = dispersion_data[:,6:end]

    obj_list = []
    for bands in nbandsbelow
        append!(obj_list, calcmincladding(freqs, bands[end], gap0lower, gap0upper))
    end
    return findmax(obj_list)
end


"""
Plot band structure. Wrapper for Brillouin's plotting function.
"""
function plotbands(calcname::String, Rs=nothing; sg=nothing, dim=nothing, mu1=nothing,
    title=nothing, band_highlights=nothing, savename=nothing, nkpoints=1000,
    hlinesy=nothing, range_fill=nothing, fillcolor="rgba(235,235,52,0.2)",
    ylim=nothing, layout=Layout(), width=nothing, height=nothing)

    if isnothing(sg)
        sg = try_parse_sgnum(calcname)    # int for space group
    end
    if isnothing(dim)
        dim = try_parse_dim(calcname)     # int for dimension
    end

    if isnothing(Rs)
        @warn "Rs was not given. Generating a random basis for plotting."
        Rs = directbasis(sg)
    end
    kp = irrfbz_path(sg, Rs)
    k_interp = Brillouin.interpolate(kp, nkpoints)

    dispersion_data = readdlm((@__DIR__)*"/output/"*calcname*"-dispersion.out", ',')
    # D = dim
    # kvecs = KVec.(eachrow(@view dispersion_data[:,2:2+(D-1)]))
    freqs = dispersion_data[:,6:end]
    nbands = size(freqs, 2)

    if !isnothing(hlinesy)
        for hliney in hlinesy
            freqs = hcat(freqs, ones((size(freqs, 1), 1)) * hliney)
        end
    end

    # band_highlights = nothing
    if !isnothing(mu1) && isnothing(band_highlights)
        band_highlights = Dict(1:mu1 => attr(color=:red, width=3))
    end
    if !isnothing(hlinesy)
        if isnothing(band_highlights)
            band_highlights = Dict(nbands+1:nbands+length(hlinesy) => attr(color=:black, width=1))
        else
            band_highlights[nbands+1:nbands+length(hlinesy)] = attr(color=:black, width=1)
        end
    end
    # println((size(k_interp), size(freqs)))

    p = Brillouin.plot(k_interp, freqs, layout, 
        title=title, band_highlights=band_highlights, ylims=ylim)

    if !isnothing(range_fill)
        local_xs = cumdists.(Brillouin.cartesianize(k_interp).kpaths)
        traces = []
        for (path_idx, local_x) in enumerate(local_xs)
            push!(traces, 
                PlotlyJS.scatter(x=vcat(local_x, reverse(local_x)),
                    y=vcat(ones(size(local_x))*range_fill[1], reverse(ones(size(local_x))*range_fill[2])),
                    fill="toself", fillcolor=fillcolor, hover="skip", xaxis="x$path_idx",
                    line=attr(color="rgba(255,255,255,0)")))
        end
        bandtrace = p.plot.data
        p = PlotlyJS.plot(append!(bandtrace, traces), p.plot.layout)
    end

    # PlotlyJS.relayout!(p, title=title)  # Doing this because older version of Brillouin doesn't take title kwarg?

    if !isnothing(savename)
        PlotlyJS.savefig(p, savename, width=width, height=height, scale=2)
    end
    return p
end

function plotbz(sgnum, Rs=nothing)
    if isnothing(Rs)
        Rs = directbasis(sgnum, Val(3)) # generate a random basis consistent with `sgnum`
    end

    cntr = centering(sgnum, 3) # will be ‘C’ for base-centered
    pRs = primitivize(Rs, cntr) # primitive direct basis
    pGs = reciprocalbasis(pRs) # primitive reciprocal basis
    cell = wignerseitz(pGs)
    kp = irrfbz_path(sgnum, Rs) # should still be the conventional basis here; think I wrote it wrong in the other email
    cartesianize!(kp) # should be the primitive reciprocal basis here

    p = Brillouin.plot(cell, kp)
    return p
end