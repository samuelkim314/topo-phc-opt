# Functions for generating and writing random PhCs
using Crystalline, MPBUtils
include("params.jl")


"""
Generate a random Fourier lattice in space group `sgnum`
"""
function init_lattice(sgnum, Dᵛ::Val{D}, cntr::Char;
                maxGs=ntuple(_->2, Val(D)), expon::Real=1.25) where D    
    # level-set surface
    flat = normscale!(modulate(levelsetlattice(sgnum, Dᵛ, maxGs)), expon)
    deleteat!(flat.orbits, 1), deleteat!(flat.orbitcoefs, 1) # trivial constant G=0 term
    pflat = primitivize(flat, cntr)

    return pflat
end

"""
Generate a random pflat, typically used as an initial point for optimization.
"""
function randpflat(sgnum::Integer, D::Integer=3;
    cntr=nothing)
    if isnothing(cntr)
        cntr = centering(sgnum, 3) # centering symbol (e.g, 'F' for face-centered, etc.)
    end

    # direct basis
    nlengthparams, nangleparams = nbasisparams(sgnum, Val(D))
    nlatticeparams = nlengthparams + nangleparams
    basisparams0 = rand(nlatticeparams)

    # generate a isosurface
    pflat0 = init_lattice(sgnum, Val(D), cntr)
    # generate some dielectric values
    filling = rand()     # Range [0, 1)
    isoval0 = MPBUtils.filling2isoval(pflat0, filling)

    # Dumb way of extracting modulation coefficients by flattening orbitcoefs
    pflat0 = Crystalline.ModulatedFourierLattice{D}([vcat(pflat0.orbits...)], [vcat(pflat0.orbitcoefs...)])    

    return isoval0, basisparams0, pflat0
end

"""
    x0, nlatticeparams, nflatparams = randinit(sgnum, D, hasinv; cntr=cntr)

Generate random PhC initialization for optimization, including filling fraction `ff`, lattice parameters, and 
Fourier lattice modulation. Also returns number of lattice parameters `nlatticeparams` and number of Fourier Lattice
parameters `nflatparams` for convenience.
"""
function randinit(sgnum::Integer, D::Integer=3, hasinv::Bool=true; cntr::Char=nothing)

    if isnothing(cntr)
        cntr = centering(sgnum, 3) # centering symbol (e.g, 'F' for face-centered, etc.)
    end

    # Count free parameters for basis lattice vectors
    nlengthparams, nangleparams = nbasisparams(sgnum, Val(D))
    nlatticeparams = nlengthparams + nangleparams
    basisparams0 = rand(nlatticeparams)

    # Generate random Fourier lattice
    pflat0 = init_lattice(sgnum, Val(D), cntr)
    # Filling factor and level-set cutoff (isoval) to generate isosurface
    ff0 = rand()    # Range [0, 1)
    isoval0 = MPBUtils.filling2isoval(pflat0, ff0)
    # Dumb way of extracting modulation coefficients by flattening orbitcoefs
    pflat0 = Crystalline.ModulatedFourierLattice{D}([vcat(pflat0.orbits...)], [vcat(pflat0.orbitcoefs...)])    

    flat0 = conventionalize(pflat0, cntr)

    modulation0 = pflat2modulation(pflat0, sgnum, cntr=cntr, hasinv=hasinv)
    nflatparams = length(modulation0)

    if hasinv
        x0 = append!([ff0], basisparams0, modulation0)
    else
        x0 = append!([ff0], basisparams0, real(modulation0), imag(modulation0))
    end

    return x0, nlatticeparams, nflatparams
end

"""
Get initial PhC lattice parameterization from an .sh file at `initpath`
"""
function sh_init(initpath, sgnum, hasinv::Bool=true; cntr::Char=nothing, abclow=0.75, abchigh=1.25)
    if isnothing(cntr)
        cntr = centering(sgnum, 3) # centering symbol (e.g, 'F' for face-centered, etc.)
    end

    pRs, pflat0, isoval0, epsin, epsout, kvecs = lattice_from_mpbparams(initpath)
    Rs = conventionalize(pRs, cntr)
    basisparams0 = basis2params(Rs, sgnum; abclow=abclow, abchigh=abchigh)
    flat0 = conventionalize(pflat0, cntr)
    ff0 = isoval2filling(flat0, isoval0)

    modulation0 = pflat2modulation(pflat0, sgnum, cntr=cntr, hasinv=hasinv)
    nlatticeparams = length(basisparams0)
    nflatparams = length(modulation0)

    if hasinv
        x0 = append!([ff0], basisparams0, modulation0)
    else
        x0 = append!([ff0], basisparams0, real(modulation0), imag(modulation0))
    end

    return x0, nlatticeparams, nflatparams
end

function sh_init_trsb(initpath, sgnum, hasinv::Bool=true; cntr::Char=nothing, abclow=0.75, abchigh=1.25)
    if isnothing(cntr)
        cntr = centering(sgnum, 3) # centering symbol (e.g, 'F' for face-centered, etc.)
    end

    pRs, pflat0, isoval0, epsin, g, epsout, kvecs = lattice_from_mpbparams_trsb(initpath)
    Rs = conventionalize(pRs, cntr)
    basisparams0 = basis2params(Rs, sgnum; abclow=abclow, abchigh=abchigh)
    flat0 = conventionalize(pflat0, cntr)
    ff0 = isoval2filling(flat0, isoval0)

    modulation0 = pflat2modulation(pflat0, sgnum, cntr=cntr, hasinv=hasinv)
    nlatticeparams = length(basisparams0)
    nflatparams = length(modulation0)

    if hasinv
        x0 = append!([ff0], basisparams0, modulation0)
    else
        x0 = append!([ff0], basisparams0, real(modulation0), imag(modulation0))
    end

    return x0, g, nlatticeparams, nflatparams
end

function write_rand_lgs(epsin=3.4757^2, 
    epsout=1.0,
    id::Integer=1,
    write_dir="input",
    sgnum::Integer=68,
    D::Integer=3,
    res::Integer=32,
    nbands::Integer=6, 
    has_tr=true)
    
    # write_dir = joinpath((@__DIR__), write_dir)
    isdir(write_dir) || mkdir(write_dir)

    cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)

    # --- find out which little groups we need to assess bandreps/symmetry vector ---
    brs  = bandreps(sgnum, D, timereversal=has_tr)
    lgs  = Crystalline.matching_littlegroups(brs, Val(D))
    plgs = primitivize.(lgs, #=modw=# false)

    # --- generate mpb-input file ---
    runtype = "all" # we don't have TE/TM in 3D; do all polarizations

    # direct basis
    Rs  = directbasis(sgnum, Val(D)) # generate a random basis consistent with `sgnum`
    pRs = primitivize(Rs, cntr)

    # generate a isosurface
    pflat = init_lattice(sgnum, Val(D), cntr)

    # generate some dielectric values
    filling = rand(0.25:0.05:0.75)

    # write to .sh file for mpb
    calcname = mpb_calcname(D, sgnum, id, res, runtype)
    pathlgs = joinpath(write_dir, calcname*".sh")
    open(pathlgs, "w") do io
        prepare_mpbcalc!(io, sgnum, pflat, pRs, filling, epsin, epsout, runtype;
                                res=res, lgs=plgs, id=id, nbands=nbands)
        # NB: if you want an mpb dispersion calculation, provide a keyword argument `kvecs`
        #     instead of `lgs` in the above. You can also supply a filename for `kvecs` if
        #     you want.
    end
    # id == 1 && println("MPB setup files written to:")
    # println("  ∙  ", filepath, ".sh")
    
    params = Dict("Rs" => Rs, "pflat" => pflat, "filling" => filling)

    return params, calcname, pathlgs
end

function write_rand_dispersion(
    epsin=3.4757^2, 
    epsout=1.0, 
    id::Integer=1, 
    write_dir="input",
    sgnum::Integer=68, 
    D::Integer=3, 
    res::Integer=32, 
    nbands::Integer=6, 
    has_tr=true)
    
    # write_dir = joinpath((@__DIR__), write_dir)
    isdir(write_dir) || mkpath(write_dir)

    cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)

    # # --- find out which little groups we need to assess bandreps/symmetry vector ---
    # brs  = bandreps(sgnum, D, timereversal=has_tr)
    # lgs  = Crystalline.matching_littlegroups(brs, Val(D))
    # plgs = primitivize.(lgs, #=modw=# false)

    # --- generate mpb-input file ---
    runtype = "all" # we don't have TE/TM in 3D; do all polarizations

    # direct basis
    Rs  = directbasis(sgnum, Val(D)) # generate a random basis consistent with `sgnum`
    pRs = primitivize(Rs, cntr)

    pflat = init_lattice(sgnum, Val(D), cntr)   # generate a isosurface

    filling = rand(0.25:0.05:0.75)  # generate some dielectric values

    kp = irrfbz_path(sgnum, Rs)
    k_interp = interpolate(kp, 1000)

    # write to .sh file for mpb
    calcname = mpb_calcname(D, sgnum, id, res, runtype)
    pathlgs = joinpath(write_dir, calcname*".sh")
    open(pathlgs, "w") do io
        prepare_mpbcalc!(io, sgnum, pflat, pRs, filling, epsin, epsout, runtype;
                                res=res, kvecs=k_interp, id=id, nbands=nbands)
        # NB: if you want an mpb dispersion calculation, provide a keyword argument `kvecs`
        #     instead of `lgs` in the above. You can also supply a filename for `kvecs` if
        #     you want.
    end
    
    params = Dict("Rs" => Rs, "pflat" => pflat, "filling" => filling)

    return params, calcname, pathlgs
end