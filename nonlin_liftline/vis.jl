using WriteVTK

# ALL BELOW IS COPIED FROM VPM CODE -----------------------------------

function mesh2Vtu(pans, verts, fname)
    cells = MeshCell[]
    for i =1:size(pans,2)
        inds = pans[:,i]
        if inds[end] == 0 # triangle
            c = MeshCell(VTKCellTypes.VTK_TRIANGLE, inds[1:end-1])
        else # quad
            c = MeshCell(VTKCellTypes.VTK_QUAD, inds)
        end
        #c = MeshCell(celltype2, inds)
        push!(cells, c)
    end
    outfile = vtk_grid(fname, verts, cells, compress=2) do vtk
    end
end

function pan2Vtu(pans, verts, panGam, panPres, panVelSurf, fname)
    cells = MeshCell[]
    for i =1:size(pans,2)
        inds = pans[:,i]
        if inds[end] == 0 # triangle
            c = MeshCell(VTKCellTypes.VTK_TRIANGLE, inds[1:end-1])
        else # quad
            c = MeshCell(VTKCellTypes.VTK_QUAD, inds)
        end
        push!(cells, c)
    end
    outfile = vtk_grid(fname, verts, cells, compress=2) do vtk
        vtk["gam"] = panGam
        vtk["pres"] = panPres
        vtk["vmag"] = sqrt.(panVelSurf[1,:].^2 .+ panVelSurf[2,:].^2 .+ panVelSurf[3,:].^2)
    end
end

function part2Vtu(partPos, fname, npart)
    if npart > 0
        x = partPos[1,1:npart]
        y = partPos[2,1:npart]
        z = partPos[3,1:npart]
    else
        npart = 1
        x = [NaN]
        y = [NaN]
        z = [NaN]
        mag = [NaN]
    end
    
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:npart]
    
    vtk_grid(fname, x, y, z, cells) do vtk
        #vtk["mag", VTKPointData()] = mag
    end
end

function line2Vtp(linePos, fname, nline)
    if nline > 0
        lines = []
        pts = zeros(3,2*nline)
        for i = 1:nline
            pts[:,(i-1)*2+1]= linePos[:,1,i]
            pts[:,(i-1)*2+2] = linePos[:,2,i]
        end
        lines = [MeshCell(PolyData.Lines(), ((i-1)*2+1, (i-1)*2+2)) for i in 1:nline]
    else
        pts = [NaN NaN; NaN NaN; NaN NaN]
        lines = [MeshCell(PolyData.Lines(), (1, 2))]
    end
    outfile = vtk_grid(fname, pts, lines) do vtk
    end
end


function norms2Vtp(norms, cpt, edgeLen, npt, fname)
    npan = size(cpt,2)
    # get average line lengths
    len = 0
    num = 0
    for i = 1:npan
        for j = 1:npt[i]
            len = len + edgeLen[j,i]
            num = num + 1
        end
    end
    len = len / num

    pts = zeros(3,2*npan)
    lines = []
    for i = 1:npan
        p1 = cpt[:,i]
        p2 = p1 .+ norms[:,i]*len
        pts[:,(i-1)*2+1]= p1
        pts[:,(i-1)*2+2] = p2
    end
    lines = [MeshCell(PolyData.Lines(), ((i-1)*2+1, (i-1)*2+2)) for i in 1:npan]
    outfile = vtk_grid(fname, pts, lines) do vtk
    end
end

function groundPlane2Vtu(planePt, planeNorm, uinf, rad, fname)
    x1 = uinf ./ norm(uinf)
    y1 = cross(planeNorm, x1)

    p1 = planePt .- x1*rad .- y1*rad
    p2 = planePt .+ x1*rad .- y1*rad
    p3 = planePt .+ x1*rad .+ y1*rad
    p4 = planePt .-x1*rad .+ y1*rad
    celltype = VTKCellTypes.VTK_QUAD
    cells = MeshCell[]
    c = MeshCell(celltype, [1;2;3;4])
    push!(cells, c)
    outfile = vtk_grid(fname, [p1 p2 p3 p4], cells, compress=2) do vtk
    end

end


