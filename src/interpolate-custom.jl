using Brillouin
import Brillouin: interpolate
using LinearAlgebra

"""
    interpolate(kp::KPath, N::Integer)
    interpolate(kp::KPath; N::Integer, density::Real) --> KPathInterpolant

Return an interpolant of `kp` with `N` points distributed approximately equidistantly
across the full **k**-path (equidistance is measured in a Cartesian metric).

Note that the interpolant may contain slightly fewer or more points than `N` (typically
fewer) in order to improve equidistance. `N` can also be provided as a keyword argument.

As an alternative to specifying the desired total number of interpolate points via `N`,
a desired density per unit (reciprocal) length can be specified via the keyword argument
`density`.
"""
function interpolate(kp::KPath{D}, N::Integer, klabskip::Symbol) where D
    kpᶜ = setting(kp) === Brillouin.CARTESIAN ? kp : Brillouin.cartesianize(kp)
    distss = map(paths(kpᶜ)) do path
        map(1:length(path)-1) do i
            norm(points(kpᶜ)[path[i]] - points(kpᶜ)[path[i+1]])
        end
    end
    totaldist = sum(dists->sum(dists), distss)

    kipaths = [Vector{Brillouin.SVector{D, Float64}}() for _ in 1:length(paths(kp))]
    labels = [Dict{Int, Symbol}() for _ in 1:length(paths(kp))]
    for (j, (path, dists)) in enumerate(zip(paths(kp), distss))
        push!(labels[j], 1 => first(path))
        for i in 1:length(path)-1
            if path[i] == klabskip || path[i+1] == klabskip
                Nᵢ = 2
            else
                Nᵢ = max(convert(Int, div(N*dists[i], totaldist, RoundUp)), 2)
            end
            append!(kipaths[j], range(points(kp)[path[i]], points(kp)[path[i+1]];
                                      length=Nᵢ)[1:end-1])
            push!(labels[j], length(kipaths[j])+1 => path[i+1])
        end
        push!(kipaths[j], points(kp)[last(path)])
    end

    return KPathInterpolant(kipaths, labels, basis(kp), Ref(setting(kp)))
end

"""
Alternative interpolation scheme to make it compatible with MPB implementation
"""
function interpolate_simple(kp::KPath{D}, N::Integer) where D
    kpᶜ = setting(kp) === Brillouin.CARTESIAN ? kp : Brillouin.cartesianize(kp)
    distss = map(paths(kpᶜ)) do path
        map(1:length(path)-1) do i
            norm(points(kpᶜ)[path[i]] - points(kpᶜ)[path[i+1]])
        end
    end

    kipaths = [Vector{Brillouin.SVector{D, Float64}}() for _ in 1:length(paths(kp))]
    labels = [Dict{Int, Symbol}() for _ in 1:length(paths(kp))]
    for (j, (path, dists)) in enumerate(zip(paths(kp), distss))
        push!(labels[j], 1 => first(path))
        for i in 1:length(path)-1
            Nᵢ = N+2
            append!(kipaths[j], range(points(kp)[path[i]], points(kp)[path[i+1]];
                                      length=Nᵢ)[1:end-1])
            push!(labels[j], length(kipaths[j])+1 => path[i+1])
        end
        push!(kipaths[j], points(kp)[last(path)])
    end

    return KPathInterpolant(kipaths, labels, basis(kp), Ref(setting(kp)))
end



# sgnum = 98
# degendim = getdegendim(sgnum)   # Degeneracy dimension
# klabskip = getdegenksym(sgnum)  # Degeneracy k-point
# cntr = centering(sgnum, D)      # centering symbol (e.g, 'F' for face-centered, etc.)
# Rs  = directbasis(sgnum, Val(D)) # generate a random basis consistent with `sgnum`
# kp = irrfbz_path(sgnum, Rs)
# k_interp = interpolate(kp, 1000, klabskip)
# # println(k_interp.labels)
# println(length(k_interp))