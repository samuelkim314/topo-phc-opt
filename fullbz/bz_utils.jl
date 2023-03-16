using Crystalline, Brillouin
using Contour, Interpolations
using LinearAlgebra # Detect I (identity matrix)
include("nodal_features_from_meshing.jl")


"""
    calc2dcentroid(freqs_cut, xs, ys)

Calculate contour around minimum of `freqs_cut` and centroid of the contour polygon

# Example
band2 = freqs[:, :, end, 2]
band3 = freqs[:, :, end, 3]
banddiff = band3 - band2
centroid = calc2dcentroid(banddiff, kxs, kys)

# Arguments
- `freqs_cut`: 2D array, typically difference between bands that we want to calculate the degeneracy
- `xs`, `ys`:   coordinates for the array `freqs_cut`
"""
function calc2dcentroid(freqs_cut, xs, ys;
        threshold=1.5)

    # Calculate contour lines
    minline = Contour.lines(Contour.contour(xs, ys, freqs_cut, threshold*minimum(freqs_cut)))
    # minline = lines(Contour.contour(kxs, kzs, freqxz, 5.0*minimum(min23)))
    contourcoordinates = Contour.coordinates(minline[1])
    contourcoordinates = [[x, y] for (x, y) in zip(contourcoordinates...)]

    cx = 0
    cy = 0
    a = 0
    npoints = size(contourcoordinates, 1)
    for j in 1:npoints-1
        xi = contourcoordinates[j][1]
        yi = contourcoordinates[j][2]
        xi1 = contourcoordinates[j%npoints+1][1]
        yi1 = contourcoordinates[j%npoints+1][2]
        aj = xi*yi1 - xi1*yi
        cx += (xi+xi1)*aj
        cy += (yi+yi1)*aj
        a += aj
    end
    a /= 2
    cx /= (6*a)
    cy /= (6*a)
    # centroid = [cx, cy, kzs[ikzweyl]]
    # println(centroid)
    return [cx, cy]
end

"""
Calculate polyhedron centroid based on its `verts` and `faces`.

Based on https://www.ma.imperial.ac.uk/~rn/centroid.pdf
"""
function calc3dcentroid(verts, faces)
    V = 0
    centroid = [0.0, 0, 0]
    for face in eachrow(faces)
        # global V, centroid
        ai, bi, ci = verts[face[1], :], verts[face[2], :], verts[face[3], :]
        nhat = cross(bi - ai, ci - ai)
        # ni = nhat ./ norm(nhat, 2)    # We don't use this variable
        V += dot(ai, nhat)
        centroid += nhat .* ((ai + bi).^2 + (bi + ci).^2 + (ci + ai).^2) / 24
    end
    V /= 6
    centroid /= (2 * V)
    return centroid
end

"""
Apply symmetry operators to Vector `point`, returns a list of coordinates
"""
function applysymmetry(point, sgnum::Integer, addinv=false)
    # Apply symmetry operators
    coords_list = [point]
    pg = pointgroup(primitivize(spacegroup(sgnum, Val(3))))
    for op in pg
        rotation(op) == I && continue
        # coords = vcat(coords, (rotation(op)'\coords')')
        newcoord = rotation(op)'\centroid
        push!(coords_list, newcoord)
    end
    coords = hcat(coords_list...)'  # Convert from list of vectors to matrix
    if addinv
        # Add inversion operator
        coords = vcat(coords, (rotation(SymOperation{3}("-x,-y,-z"))'\coords')')
    end
    return coords
end


"""
    Interpolate bands at point `coord` to get the Weyl point frequency
"""
function getweylfreq(coord, freqs, kxs, kys, kzs, nbandsbelow)
    itpbelow = LinearInterpolation((kxs, kys, kzs), freqs[:, :, :, nbandsbelow])
    itpabove = LinearInterpolation((kxs, kys, kzs), freqs[:, :, :, nbandsbelow+1])
    weylfreqbelow = itpbelow(coord...)
    weylfreqabove = itpabove(coord...)
    weylfreq = (weylfreqbelow + weylfreqabove) / 2
    println((weylfreqbelow, weylfreqabove, weylfreq, (weylfreqabove - weylfreqbelow)/weylfreq))
    return weylfreq
end


function isosurfacemesh(data, isoval, nbandsbelow, cartesianize=true, addinv=true;
    invertfaces=true)
    # Plot Fermi pockets
    freqs, (kxs, kys, kzs), (Nkx, Nky, Nkz), (Δx, Δy, Δz) = unpack_mpb_data(data)
    origin=SVector(minimum(kxs), minimum(kys), minimum(kzs))
    widths=SVector(Δx, Δy, Δz)
    # verts, faces = Meshing.isosurface(freqs4, MarchingTetrahedra(iso=isoval, eps=1e-5), origin=origin, widths=widths)
    # Gstemp=ReciprocalBasis{3}([1.,0,0], [0,1.,0], [0,0,1.])
    # verts, faces = _mesh_to_cartesian(verts, faces, Gstemp)
    freqs0 = reshape(freqs[:,nbandsbelow], Nkx, Nky, Nkz)
    freqs1 = reshape(freqs[:,nbandsbelow+1], Nkx, Nky, Nkz)
    verts0, faces0 = nodal_features_as_mesh_symmetry_extended(freqs0, isoval, sgnum; origin=origin, widths=widths, addinv=addinv)
    verts1, faces1 = nodal_features_as_mesh_symmetry_extended(freqs1, isoval, sgnum; origin=origin, widths=widths, addinv=addinv)
    verts = vcat(verts0, verts1)
    faces = vcat(faces0, faces1.+size(verts0, 1))

    if cartesianize
        # Convert vertices to Cartesian basis
        for (i, vert) in enumerate(eachrow(verts))
            verts[i, :] = Brillouin.cartesianize(vert, Gs)
        end
    end

    if invertfaces
        faces[:, 2], faces[:, 3] = faces[:, 3], faces[:, 2]
    end

    return verts, faces
end
