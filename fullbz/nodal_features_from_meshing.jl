using StaticArrays, Meshing
using PyPlot # need to load since _mesh_to_cartesian is conditional on PyPlot being load
using Crystalline: _mesh_to_cartesian, ReciprocalBasis
using PlotlyJS
using LinearAlgebra: norm
include("brillouin_zones.jl"); using Main.BrillouinZone


function _isosurface(df, level, 
    Algo::UnionAll=MarchingTetrahedra;
    Gs::ReciprocalBasis{3}=ReciprocalBasis{3}([1.,0,0], [0,1.,0], [0,0,1.]),
    eps=1e-5, origin=SVector(0.0,0.0,0.0), widths=SVector(1.0,1.0,1.0))

    verts, faces = Meshing.isosurface(df, Algo(iso=level, eps=eps), origin=origin, widths=widths)
    verts′, faces′ = _mesh_to_cartesian(verts, faces, Gs)
end

function unpack_mpb_data(data)

    freqs = sort!(data[:,6:end], dims=2)
    kvs = collect(eachrow(data[:,2:4]))
    kxs, kys, kzs = sort(unique(getindex.(kvs, 1))), sort(unique(getindex.(kvs, 2))), sort(unique(getindex.(kvs, 3)))
    Δx, Δy, Δz = -(-)(extrema(kxs)...), -(-)(extrema(kys)...), -(-)(extrema(kzs)...)
    Nkx, Nky, Nkz = length(kxs), length(kys), length(kzs)

    return freqs, (kxs, kys, kzs), (Nkx, Nky, Nkz), (Δx, Δy, Δz)
end

function nodal_features_as_mesh(data;
            bands = (2,3), grace_factor = 15.0, level = nothing,
            Gs::ReciprocalBasis{3}=ReciprocalBasis{3}([1.,0,0], [0,1.,0], [0,0,1.]))
    
    # unpacking data
    freqs, (kxs, kys, kzs), (Nkx, Nky, Nkz), (Δx, Δy, Δz) = unpack_mpb_data(data)

    # finding nodal lines
    df = reshape(abs.(freqs[:,bands[2]] .- freqs[:,bands[1]]), (Nkx, Nky, Nkz))
    if level === nothing && grace_factor !== nothing
        δ = minimum(df)
        println(findmin(df))
        level = grace_factor*δ
    end
    println( (Δx, Δy, Δz))
    verts, faces = _isosurface(df, level; Gs=Gs, 
                            origin=SVector(minimum(kxs), minimum(kys), minimum(kzs)),
                            widths=SVector(Δx, Δy, Δz) )

    return verts, faces, freqs, df, (kxs, kys, kzs), (Nkx, Nky, Nkz)
end


function nodal_features_as_mesh_symmetry_extended(df, level, sgnum;
            Gs::Union{ReciprocalBasis{3}, Nothing}=nothing, 
            algorithm::UnionAll=MarchingTetrahedra,
            origin=SVector(0.0,0.0,0.0), widths=SVector(1.0,1.0,1.0),
            addinv=false)

    verts⁰, faces⁰ = _isosurface(df, level, algorithm, origin=origin, widths=widths)

    pg = pointgroup(primitivize(spacegroup(sgnum, Val(3))))
    verts, faces = copy(verts⁰), copy(faces⁰)
    # use symmetries to extend to full BZ
    for op in pg
        # println((op, rotation(op), Crystalline.rotation_order(op)))
        rotation(op) == I && continue
        facesop = faces⁰ .+ size(verts,1)
        if Crystalline.rotation_order(op) == -2
            facesop[:, 2], facesop[:, 3] = facesop[:, 3], facesop[:, 2]
        end
        faces = vcat(faces, facesop)
        verts = vcat(verts, (rotation(op)'\verts⁰')')
    end
    if addinv && !isempty(faces) 
        # verts0, faces0 = copy(verts), copy(faces)
        facesinv = faces .+ size(verts, 1)
        facesinv[:, 2], facesinv[:, 3] = facesinv[:, 3], facesinv[:, 2]
        faces = vcat(faces, facesinv)
        verts = vcat(verts, (rotation(SymOperation{3}("-x,-y,-z"))'\verts')')
    end

    # TODO: ought to clean up trivial repetitions in verts afterwards... bit tedious

    #=
    TODO: Enable and add kwarg access point?
    # enforce periodicity
    for i in -1:1, j in -1:1
    i==0 && j==0 && continue
    global verts, faces, verts⁰, faces⁰
    faces = vcat(faces, faces⁰ .+ size(verts,1))
    verts = vcat(verts, verts⁰ .+ [i j 0.0])
    end
    # remove verts/face outside original "inner" BZ
    vertsvec, facesvec = collect(eachrow(verts)), collect(eachrow(faces))
    cutval = 0.5
    killfaces_idxs = findall(idxs->any(i->any(x->abs(x)>cutval, vertsvec[i]), idxs), facesvec)
    deleteat!(facesvec, killfaces_idxs)
    faces = hcat(facesvec...)'
    =#

    if Gs !== nothing    # convert to cartesian coordinates
        Gm = hcat(Gs...)
        verts = (Gm*verts')' 
    end

    return verts, faces
end

"""
    plot_nodal_lines_in_bz(verts, faces, sgnum::Integer, Rs::DirectBasis{3})

## Output

`traces` and `layout`. Plot these objects via

```
p = PlotlyJS.Plot(traces, layout) # create the underlying plot object (::Plot)
P = PlotlyJS.plot(p)              # visualize the plot object in a ::SyncPlot
PlotlyJS.display_blink(P)         # sometimes the visualizion fails to appear: force it to
```
We return `traces` and `layout` (rather than just plotting them directly) in case the callee
wants to add additional traces (i.e. graphs) before executing the plot.

## Input

- `Rs`: *primitive* direct basis (ITA conventions)
"""
function plot_nodal_lines_in_bz(verts, faces, sgnum::Integer, Rs::DirectBasis{3})
    # Nodal lines/mesh
    tnl = mesh3d(x=verts[:,1], y=verts[:,2], z=verts[:,3], 
                i=faces[:,1].-1, j=faces[:,2].-1, k=faces[:,3].-1, 
                color="rgb(33,64,154)", hovertext="nodal (out)line")

    # BZ
    Gs = reciprocalbasis(Rs); Gm = hcat(Gs...)
    sgnum ≠ 0 ? (cRs = conventionalize(Rs, centering(sgnum, 3))) : (cRs = Rs) # let sgnum=0 be sentinel value for skipping this and doing paralleliped
    lab2kv, cncts = BrillouinZone.brillouin_zone(sgnum, cRs)
    # ... lines
    tbz = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, length(cncts))
    for (idx, cnct) in enumerate(cncts)
        kvpath = [Gm*lab2kv[lab] for lab in cnct]
        tbz[idx] = PlotlyJS.scatter3d(
            x=getindex.(kvpath, 1), y=getindex.(kvpath, 2), z=getindex.(kvpath, 3),
            line=attr(color="rgb(0,0,0)", width=2), mode="lines",
            hovertext="BZ", hoverinfo="text"
            )
    end
    # ... "origo"
    tor = PlotlyJS.scatter3d(x=[0], y=[0], z=[0], mode="markers",
        marker=attr(color="rgb(0,0,0)", size=5, symbol="cross"), 
        hovertext="(0,0,0)", hoverinfo="text")

    # Reciprocal lattice vectors
    tgs    = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 3)
    tgtips = Vector{GenericTrace{Dict{Symbol,Any}}}(undef, 3)
    gcol   = "rgb(96,96,96)"
    startfac = .525
    extG   = .25*sum(norm.(Gs))/3
    for (idx, G) in enumerate(Gs)
        G′ = G./norm(G)
        G₀ = G*startfac
        G₁ = G₀ .+ G′.*extG
        name = "<b>G</b><sub>$(idx)</sub>"
        tgs[idx] = PlotlyJS.scatter3d(
            x=[G₀[1],G₁[1]], y=[G₀[2],G₁[2]], z=[G₀[3],G₁[3]];
            mode="lines", line=attr(color=gcol, width=4.5), 
            hovertext=name, hoverinfo="text")
        tgtips[idx] = PlotlyJS.cone(
            x=[G₁[1]], y=[G₁[2]], z=[G₁[3]], u=[G′[1]], v=[G′[2]], w=[G′[3]],
            sizeref=.5, showscale=false, anchor="tail", colorscale=[[0, gcol], [1, gcol]],
            hovertext=name, hoverinfo="text")
    end

    # Settings
    layout=Layout(
        showlegend=false,
        scene=attr(
            xaxis=attr(tickvals=[], zeroline=false,
                    showgrid=false, showbackground=false,
                    title=attr(text="")#"k<sub>x</sub>")
                    ),
            yaxis=attr(tickvals=[], zeroline=false,
                    showgrid=false, showbackground=false,
                    title=attr(text="")#"k<sub>y</sub>")
                    ),
            zaxis=attr(tickvals=[], zeroline=false,
                    showgrid=false, showbackground=false,
                    title=attr(text=""),#"k<sub>z</sub>")   
                    ),
            aspectmode = "data",
            ),
        margin=attr(l=0, r=0, b=0, t=0),
        autosize=false,            
    )

    traces = [tnl, tbz..., tor, tgs..., tgtips...]

    return traces, layout
end