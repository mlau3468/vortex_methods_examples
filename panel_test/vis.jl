
function panels2vtk(panel_list, rr, fname)
    celltype = VTKCellTypes.VTK_QUAD
    cells = MeshCell[]
    cdata = Float32[]
    pts = rr
    for i =1:size(panel_list,1)
        inds = panel_list[i].ee
        c = MeshCell(celltype, inds)
        push!(cells, c)
    end
    outfile = vtk_grid(fname, pts, cells, compress=2) do vtk
    end
end

function eerr2vtk(ee, rr, fname)
    celltype = VTKCellTypes.VTK_QUAD
    cells = MeshCell[]
    cdata = Float32[]
    pts = rr
    for i =1:size(ee,2)
        inds = ee[:,i]
        c = MeshCell(celltype, inds)
        push!(cells, c)
    end
    outfile = vtk_grid(fname, pts, cells, compress=2) do vtk
    end
end
