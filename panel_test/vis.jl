
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

npoints = 100
x = rand(npoints)
y = rand(npoints)
z = rand(npoints)
pressure = rand(npoints)
temp = rand(npoints)
cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:npoints]
vtk_grid("./points", x, y, z, cells) do vtk
    vtk["pressure", VTKPointData()] = pressure
    vtk["temp", VTKPointData()] = temp
end

npoints = size(particles_list, 1)
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
