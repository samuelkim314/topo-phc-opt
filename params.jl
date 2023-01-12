#=
Functions for processing photonic crystal parameters for optimization.
=#

using Crystalline
using Crystalline: AbstractFourierLattice
using MPBUtils

°(φ::Real) = deg2rad(φ)


"""
Get the number of free parameters for the basis lattice vectors for a given space group
"""
function nbasisparams(sgnum::Integer, Dᵛ::Val{D}=Val(3);) where D
    system = crystalsystem(sgnum, D)
    if D == 1
        return 0, 0
    elseif D == 2
        if     system == "square"      # a=b & γ=90° (free: a)
            return 0, 0
        elseif system == "rectangular" # γ=90° (free: a,b)
            return 1, 0         
        elseif system == "hexagonal"   # a=b & γ=120° (free: a)
            return 0, 0
        elseif system == "oblique"     # no conditions (free: a,b,γ)
            return 1, 1
        else 
            throw(DomainError(system))
        end

    elseif D == 3
        if     system == "cubic"        # a=b=c & α=β=γ=90° (free: a)
            return 0, 0
        elseif system == "hexagonal" || # a=b & α=β=90° & γ=120° (free: a,c)
               system == "trigonal"    
            return 1, 0
        elseif system == "tetragonal"   # a=b & α=β=γ=90° (free: a,c)
            return 1, 0
        elseif system == "orthorhombic" # α=β=γ=90° (free: a,b,c)
            return 2, 0
        elseif system == "monoclinic"   # α=γ=90° (free: a,b,c,β≥90°)
            return 2, 1
        elseif system == "triclinic"    # no conditions (free: a,b,c,α,β,γ)
            return 2, 3
        else 
            throw(DomainError(system))
        end        
    else 
        _throw_invaliddim(D)
    end
end


"""
Generate basis lattice vectors from the free parameters.

Note that if the lattice is cubic, params will be ignored.
"""
function params2basis(params::Vector, sgnum::Integer, Dᵛ::Val{D}=Val(3);
            abclow=0.5, abchigh=2.0) where D
    system = crystalsystem(sgnum, D)
    αβγlow = °(30)
    αβγhigh = °(150)
    if D == 1
        return 0, 0
    elseif D == 2
        throw("2D basis not implemented")

    elseif D == 3
        if system == "cubic"        # a=b=c & α=β=γ=90° (free: a)
            a = b = c = 1.0
            α = β = γ = °(90)
        elseif system == "hexagonal" || # a=b & α=β=90° & γ=120° (free: a,c)
               system == "trigonal"    
            a = b = 1.0;        
            c = params[1] * (abchigh - abclow) + abclow
            α = β = °(90);      γ = °(120)
        elseif system == "tetragonal"   # a=b & α=β=γ=90° (free: a,c)
            a = b = 1.0;        
            c = params[1] * (abchigh - abclow) + abclow
            α = β = γ = °(90)
        elseif system == "orthorhombic" # α=β=γ=90° (free: a,b,c)
            a = 1.0;
            b = params[1] * (abchigh - abclow) + abclow
            c = params[2] * (abchigh - abclow) + abclow
            α = β = γ = °(90)
        elseif system == "monoclinic"   # α=γ=90° (free: a,b,c,β≥90°)
            a = 1.0;
            b = params[1] * (abchigh - abclow) + abclow
            c = params[2] * (abchigh - abclow) + abclow
            α = γ = °(90);      
            β = params[3] * (αβγhigh - °(90)) + °(90)
        elseif system == "triclinic"    # no conditions (free: a,b,c,α,β,γ)
            a = 1.0;
            U = _Uniform(αβγlims...)
            α, β, γ = rand(U), rand(U), rand(U)
            throw("Triclinic basis not implemented")
        else 
            throw(DomainError(system))
        end        
        return crystal(a,b,c,α,β,γ)

    else 
        _throw_invaliddim(D)
    end
end

"""
Generate vector of free parameters from the basis lattice vectors.
"""
function basis2params(Rs, sgnum::Integer, Dᵛ::Val{D}=Val(3);
            abclow=0.5, abchigh=2.0) where D
    # Make sure to pass in conventionalized basis!
    αβγlow = °(30)
    αβγhigh = °(150)
    system = crystalsystem(sgnum, D)
    if D == 1
        return 0, 0
    elseif D == 2
        throw("2D basis not implemented")

    elseif D == 3
        if system == "cubic"        # a=b=c & α=β=γ=90° (free: a)
            return nothing
        elseif system == "hexagonal" || # a=b & α=β=90° & γ=120° (free: a,c)
               system == "trigonal"    
            return [(Rs[3][3] - abclow) / (abchigh - abclow)]
        elseif system == "tetragonal"   # a=b & α=β=γ=90° (free: a,c)
            return [(Rs[3][3] - abclow) / (abchigh - abclow)]
        elseif system == "orthorhombic" # α=β=γ=90° (free: a,b,c)
            p2 = (Rs[2][2] - abclow) / (abchigh - abclow)
            p3 = (Rs[3][3] - abclow) / (abchigh - abclow)
            return p2, p3
        elseif system == "monoclinic"   # α=γ=90° (free: a,b,c,β≥90°)
            R1 = Rs[1]
            R3 = Rs[3]
            β = acos(dot(R1, R3) / (norm(R1) * norm(R3)))
            p2 = (Rs[2][2] - abclow) / (abchigh - abclow)
            p3 = (norm(R3) - abclow) / (abchigh - abclow)
            β = (β - °(90)) / (αβγhigh - °(90))
            return p2, p3, β
        elseif system == "triclinic"    # no conditions (free: a,b,c,α,β,γ)
            throw("Triclinic basis not implemented")
        else 
            throw(DomainError(system))
        end        

    else 
        _throw_invaliddim(D)
    end
end

"""
Convert NLopt parameter into description of photonic crystal that can be used to run MPB.

if isnothing(basisparams), then assumes that basis params are part of x.

x should be structured [isoval; basisparams; flatparams] or [ff; basisparams; flatparams]
"""
function params2flat(x::Vector, sgnum::Integer; 
    cntr::Char=nothing, hasinv::Bool=nothing, abclow=0.5, abchigh=2.0, basisparams=nothing)

    if isnothing(cntr)
        cntr = centering(sgnum, 3) # centering symbol (e.g, 'F' for face-centered, etc.)
    end
    if isnothing(hasinv) 
        hasinv = SymOperation{3}("-x,-y,-z") in spacegroup(sgnum, D).operations
    end

    if isnothing(basisparams)
        nlengthparams, nangleparams = nbasisparams(sgnum, Val(3))
        nlatticeparams = nlengthparams + nangleparams
    else
        nlatticeparams = 0
    end
    
    # Separate out isoval, basis, and modulation parameters
    isoval = x[1]
    if isnothing(basisparams)
        # This still works even if nlatticeparams=0 (which implies cubic lattice)
        basisparams = x[2:2+nlatticeparams-1]
    end
    Rs = params2basis(basisparams, sgnum, abclow=abclow, abchigh=abchigh)
    pRs = primitivize(Rs, cntr)

    if hasinv
        modulation = x[2+nlatticeparams:end]
    else
        # First half of parameters (modulation_split) are the real parts, second half
        # are imaginary parts of modulation coefficients
        modulation_split = x[2+nlatticeparams:end]
        nmod = length(modulation_split) ÷ 2
        modulation = modulation_split[1:nmod] + modulation_split[nmod+1:end]*im
    end

    # Fourier lattice
    uflat = levelsetlattice(sgnum, 3)
    deleteat!(uflat.orbits, 1), deleteat!(uflat.orbitcoefs, 1) # trivial constant G=0 term
    uflat = primitivize(uflat, cntr)
    modulated_orbitcoefs = uflat.orbitcoefs.*modulation
    flat = Crystalline.ModulatedFourierLattice{3}(uflat.orbits, modulated_orbitcoefs)

    return isoval, pRs, flat
end

"""
Convert flat into a modulation vector, the coefficients of the orbits

Note that `pflat` is a flattened version, which can be done by
```pflat0 = Crystalline.ModulatedFourierLattice{D}([vcat(pflat0.orbits...)], [vcat(pflat0.orbitcoefs...)])```
TODO: Make this more robust
"""
function pflat2modulation(pflat::AbstractFourierLattice, sgnum::Integer;
    cntr::Char=nothing, hasinv::Bool=nothing)
    if isnothing(cntr)
        cntr = centering(sgnum, 3) # centering symbol (e.g, 'F' for face-centered, etc.)
    end
    if isnothing(hasinv) 
        hasinv = SymOperation{3}("-x,-y,-z") in spacegroup(sgnum, D).operations
    end
    # Primitive unity lattice - defines the structure and relationships of orbits
    uflat = levelsetlattice(sgnum, 3)
    deleteat!(uflat.orbits, 1), deleteat!(uflat.orbitcoefs, 1) # trivial constant G=0 term
    uflat = primitivize(uflat, cntr)

    # Extract the modulation coefficients from the initial condition
    # If it has inversion symmetry, we can constrain modulation coefficients to be real
    if hasinv
        modulation = zeros(Float64, length(uflat.orbitcoefs))
    else
        modulation = zeros(ComplexF64, length(uflat.orbitcoefs))
    end
    for i=1:length(uflat.orbitcoefs)
        orbit_i = uflat.orbits[i][1]
        ucoef_i = uflat.orbitcoefs[i][1]
        flat_ind = findfirst(x -> x == orbit_i, pflat.orbits[1])
        if hasinv
            modulation[i] = real(pflat.orbitcoefs[1][flat_ind]) / real(ucoef_i)
        else
            modulation[i] = (pflat.orbitcoefs[1][flat_ind] / ucoef_i)
        end
    end

    # Restructuring to group connected orbits - not necessary unless we want to apply (un)normscale
    # pflat0 = Crystalline.ModulatedFourierLattice{D}(uflat.orbits, uflat.orbitcoefs.*modulation0)

    return modulation
end

"""
Generate a random Fourier lattice given `sgnum`
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