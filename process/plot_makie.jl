# Plotting utilities using GLMakie
using Crystalline
using PyPlot    # This is required for mesh_3d_levelsetlattice
using GLMakie
using Crystalline: AbstractFourierLattice

function plot_pc(flat::Crystalline.AbstractFourierLattice, Rs::DirectBasis{D}, 
        isoval::Union{Real, Nothing}=nothing, N::Integer=(D==2 ? 100 : 20),
        isocaps::Bool=false) where D

    xyz = range(-.5, .5, length=N)
    vals = Crystalline.calcfouriergridded(xyz, flat, N)
    if isocaps
        vals[1,:,:] .= 0
        vals[:,1,:] .= 0
        vals[:,:,1] .= 0
        vals[end,:,:] .= 0
        vals[:,end,:] .= 0
        vals[:,:,end] .= 0
    end
    verts, faces = Crystalline.mesh_3d_levelsetlattice(vals, isoval, Rs)

    f = GLMakie.mesh(verts, faces, color=verts[:, 1])
    # GLMakie.mesh(X, T, color=X[:, 2], shading=true, figure=(resolution=(1000, 1000),))
    # hidedecorations!(ax)
    # hidespines!(ax)

    return f
end

function plot_pc_rotated(flat::Crystalline.AbstractFourierLattice, Rs::DirectBasis{D}, 
        isoval::Union{Real, Nothing}=nothing, N::Integer=(D==2 ? 100 : 20)) where D

    xyz = range(-.5, .5, length=N)
    vals = Crystalline.calcfouriergridded(xyz, flat, N)
    verts, faces = Crystalline.mesh_3d_levelsetlattice(vals, isoval, Rs)

    fig = GLMakie.Figure()
    ax = Axis3(fig[1, 1], azimuth=-0.2pi, aspect=:data, perspectiveness=0.2)
    # ax = Axis3(fig[1, 1])
    hidedecorations!(ax)
    hidespines!(ax)
    f = GLMakie.mesh!(ax, verts, faces, color=verts[:, 1])
    return fig
end

function plot_bz(sgnum::Integer)
    # n = 39
    kp = irrfbz_path(sgnum, directbasis(sgnum,3))
    c = wignerseitz(basis(kp))
    # f,ax,p = plot(c)
    # plot!(kp)
    # f
end