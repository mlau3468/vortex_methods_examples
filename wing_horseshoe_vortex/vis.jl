using WriteVTK
"""
    vis_mesh(pan_con, pan_vert, npan, fname)
Visualize mesh defined by node connectivity pan_con and list of vertices pan_vert using VTK.
"""
function vis_mesh(pan_con, pan_vert, npan, fname)
    cells = MeshCell[]
    if npan > 0
        for i =1:npan
            inds = pan_con[:,i]
            if inds[end] == 0 # triangle
                c = MeshCell(VTKCellTypes.VTK_TRIANGLE, inds[1:end-1])
            else # quad
                c = MeshCell(VTKCellTypes.VTK_QUAD, inds)
            end
            push!(cells, c)
        end
        outfile = vtk_grid(fname, pan_vert, cells, compress=2) do vtk
        end
    else
        verts = NaN .* ones(3,4)
        inds = [1;2;3;4]
        c = MeshCell(VTKCellTypes.VTK_QUAD, inds)
        push!(cells, c)
        outfile = vtk_grid(fname, pan_vert, cells, compress=2) do vtk
        end
    end
end

function vis_liftline(pan_con, pan_vert, fname)
    cells = MeshCell[]
    npan = size(pan_con, 2)
    for i =1:npan
        # quad on wing
        c = MeshCell(VTKCellTypes.VTK_QUAD, pan_con[[2;3;4;5],i])
        push!(cells, c)
        # quad on wake
        c = MeshCell(VTKCellTypes.VTK_QUAD, pan_con[[1;2;5;6],i])
        push!(cells, c)
    end
    outfile = vtk_grid(fname, pan_vert, cells, compress=2) do vtk
    end
end

function vis_normcpt(pan_cpt, pan_norm, fname)
    x = pan_cpt[1,:]
    y = pan_cpt[2,:]
    z = pan_cpt[3,:]
    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:size(pan_cpt,2)]
    vtk_grid(fname, x, y, z, cells) do vtk
        vtk["norm",  VTKPointData()] = pan_norm
    end
end