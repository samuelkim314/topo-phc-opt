using Crystalline
using Crystalline: formatirreplabel
using SymmetryBases
using PhotonicBandConnectivity
using MPBUtils
include("topology.jl")

# temporary
using Crystalline: label, formatirreplabel, symvec2string

"""
Calculate band summaries with topological properties.
See https://github.com/thchr/MPBUtils.jl for more information.
"""
function calcbandsummaries(sgnum, calcname::String; dir::String="output/", timereversal=true)
    D = 3
    symeigsd, lgd = read_symdata(calcname; sgnum=sgnum, D=D, dir=dir)
    # --- fix singular photonic symmetry content at Γ, ω=0 ---
    lgs = littlegroups(sgnum, Val(D)) 
    fixup_gamma_symmetry!(symeigsd, lgs["Γ"])
    lgirsd = pick_lgirreps(lgd; timereversal=timereversal)
    brs = bandreps(sgnum; timereversal=timereversal)
    # --- analyze connectivity and topology of symmetry data ---
    summaries = analyze_symmetry_data(symeigsd, lgirsd, brs)
    return summaries
end

"""
Calculate band topology in a simpler form than `calcbandsummaries`. Returns a dict with `UnitRange` corresponding to the
band grouping as a key, and the topological index as the value.
"""
function calcbandtopo(sgnum, calcname::String;
        dir::String="output/")
    summaries=calcbandsummaries(sgnum, calcname; dir=dir)
    dict_topology = Dict()
    for summary in summaries
        dict_topology[summary.bands] = summary.topology
    end
    return dict_topology
end

"""
Calculate band topology, grouping the bottom bands. I.e. if bands 1:3 and 4:6 are each connected, then this will
calculate the topology for 1:6.
"""
function checknontriviallowest(sgnum, calcname::String;
    dir::String="output/")

    summaries=calcbandsummaries(sgnum, calcname; dir=dir)
    dict_topology = Dict()

    for i in 1:length(summaries)
        sumbands = sum(summaries[1:i])
        if !all(x->x==0, sumbands.indicators)   # Not trivial band grouping
            dict_topology[sumbands.bands] = sumbands.indicators
        end
    end

    return dict_topology
end

function checkdegen(sgnum, calcname)
    # --- setup info ---
    has_tr = true

    # --- hilbert/ebr bases ---
    sb, brs = compatibility_basis(sgnum, 3; timereversal=has_tr)

    # --- load and process data ---
    bandirsd, lgirsd = extract_individual_multiplicities(
                            calcname,
                            timereversal=has_tr,
                            dir = "output/",
                            atol=2e-2,
                            D=3)
    length(lgirsd) ≠ length(sb.klabs) && error("missing k-point data")

    # CharacterTable(get_lgirreps(222, 3)["R"])
    # println(bandirsd["R"])    # Print degeneracies at point R

    irlabs = Dict(klab => formatirreplabel.(label.(lgirs)) for (klab, lgirs) in lgirsd)
    prettybands = Dict(klab => [bands => symvec2string(n, irlabs[klab]; braces=false)
                     for (bands, n) in bandirs]         for (klab, bandirs) in bandirsd)
    # println(irlabs)

    degenbands = UnitRange[]
    for (band, char) in prettybands[getdegenklab(sgnum)]
        if length(band) == getdegendim(sgnum) && char == getdegenchar(sgnum)
            push!(degenbands, band)
        end
    end
   
    return degenbands
end

"""
Check symmetry eigenvalues to see if the bands allow for a hinge state. Only returns lowest set of bands.
    
Arguments
    - sgnum: space group number
    - calcname: file path without "output/" prefix
    - topoindex: topological index to check for
"""
function checkhinge(sgnum::Integer, calcname::String; 
        topoindex=[0, 6], dir::String="output/")
    timereversal = false
    D = 3
    symeigsd, lgd = read_symdata(calcname; sgnum=sgnum, D=D, dir=dir)
    # --- fix singular photonic symmetry content at Γ, ω=0 ---
    fixup_gamma_symmetry!(symeigsd, lgd["Γ"])
    lgirsd = pick_lgirreps(lgd; timereversal=timereversal)
    brs = bandreps(sgnum; timereversal=timereversal)
    # --- analyze connectivity and topology of symmetry data ---
    summaries = analyze_symmetry_data(symeigsd, lgirsd, brs)

    for i in 1:length(summaries)
        sumbands = sum(summaries[1:i])
        if sumbands.indicators == topoindex
            return sumbands.bands
        end
    end
    return nothing
end

"""
Check symmetry eigenvalues to see if the bottom bands allow for a particular topological index
    
Arguments
    - sgnum: space group number
    - calcname: file path without "output/" prefix
    - topoindex: topological index to check for
"""
function checktopoindex(sgnum::Integer, calcname::String; 
        topoindex=[0, 6], dir::String="output/")
    timereversal = false
    D = 3
    symeigsd, lgd = read_symdata(calcname; sgnum=sgnum, D=D, dir=dir)
    # --- fix singular photonic symmetry content at Γ, ω=0 ---
    fixup_gamma_symmetry!(symeigsd, lgd["Γ"])
    lgirsd = pick_lgirreps(lgd; timereversal=timereversal)
    brs = bandreps(sgnum; timereversal=timereversal)
    # --- analyze connectivity and topology of symmetry data ---
    summaries = analyze_symmetry_data(symeigsd, lgirsd, brs)

    bands_list = []
    for i in 1:length(summaries)
        sumbands = sum(summaries[1:i])
        if sumbands.indicators == topoindex
            push!(bands_list, sumbands.bands)
        end
    end
    return bands_list
end

"""
Check symmetry eigenvalues to see if the bottom bands allow for a trivial gap
"""
function checktriviallowest(sgnum, calcname; Nmax=8, dir::String="output/")
    # --- setup info ---
    has_tr = false

    # --- hilbert/ebr bases ---
    sb, brs = compatibility_basis(sgnum, 3; timereversal=has_tr)
    B = matrix(brs, true)
    F = smith(B)

    # TODO: Replace with read_symdata to get symeigsd
    # And then analyze_symmetry_data to get band info

    # --- load and process data ---
    bandirsd, lgirsd = extract_individual_multiplicities(
                            calcname,
                            timereversal=has_tr,
                            dir=dir,
                            atol=2e-2,
                            D=3)
    length(lgirsd) ≠ length(sb.klabs) && error("missing k-point data")

    # extract the _potentially_ separable symmetry vectors `ns` and their band-ranges `bands`
    bands = nothing
    ns = nothing
    try
        bands, ns = extract_candidate_symmetryvectors(bandirsd, lgirsd, brs)
    catch e
        if e == ErrorException("found no isolable band candidates")
            @warn "found no isolable band candidates"
            return nothing
        else
            throw(e)
        end
    end

    if isempty(ns)
        @warn "found no isolable band candidates"
        return []
    end

    bands_trivial = UnitRange{Integer}[]

    for N in 1:Nmax
        # find symmetry vector `n′` of band grouping `band′`
        n′ = zero(first(ns))
        lower_band_seen = false
        band′ = 1:N
        band′_check = 0:0
        for (band, n) in zip(bands, ns)
            if !lower_band_seen
                if minimum(band) != minimum(band′)
                    continue
                else
                    band′_check = band
                    lower_band_seen = true
                end
            else
                if maximum(band) > maximum(band′)
                    break
                end
            end
            n′ = n′ + n
            band′_check = first(band′_check):last(band)
        end

        if band′_check == band′ && is_transverse_bandstruct(n′, sb, lgirsd["Γ"], F)
            ν,_ = indicators_singular(n′, sb, lgirsd["Γ"], B)
            if all(x->x==0, ν)
                # We have a winner
                push!(bands_trivial, band′)
            end
        end
    end

    return bands_trivial
end

