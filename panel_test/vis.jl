
function panels2vtk(panel_list, rr, fname, vis_dir)
    celltype = VTKCellTypes.VTK_QUAD
    cells = MeshCell[]
    cdata = Float32[]
    pts = rr
    mag = zeros(size(panel_list,1))
    for i =1:size(panel_list,1)
        inds = panel_list[i].ee
        c = MeshCell(celltype, inds)
        push!(cells, c)
        mag[i] = panel_list[i].mag[1]
    end
    outfile = vtk_grid(vis_dir * fname, pts, cells, compress=2) do vtk
        vtk["mag"] = mag
    end
end

function eerr2vtk(ee, rr, fname, vis_dir)
    celltype = VTKCellTypes.VTK_QUAD
    cells = MeshCell[]
    cdata = Float32[]
    pts = rr
    for i =1:size(ee,2)
        inds = ee[:,i]
        c = MeshCell(celltype, inds)
        push!(cells, c)
    end
    outfile = vtk_grid(vis_dir * fname, pts, cells, compress=2) do vtk
    end
end

function particles2vtk(particles_list, fname, vis_dir)

npoints = size(particles_list, 1)
if npoints > 0
x = zeros(size(particles_list, 1))
y = zeros(size(particles_list, 1))
z = zeros(size(particles_list, 1))
mag = zeros(size(particles_list, 1))
for i = 1:size(particles_list, 1)
    x[i] = particles_list[i].center[1]
    y[i] = particles_list[i].center[2]
    z[i] = particles_list[i].center[3]
    mag[i] = particles_list[i].mag[1]
end

cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:npoints]

vtk_grid(vis_dir * fname, x, y, z, cells) do vtk
    vtk["mag", VTKPointData()] = mag
end
end
end


function debugMatrices(A, RHS, sol, it, debug_dir)
    it = Int(it)
    # LU decomposition of A to check with DUST
    A_LU = factorize2(A)

    writedlm(debug_dir*"A_$it.csv", A)
    writedlm(debug_dir*"A_LU_$it.csv", A_LU)
    writedlm(debug_dir*"B_$it.csv", RHS)
    writedlm(debug_dir*"res_$it.csv", sol)
end