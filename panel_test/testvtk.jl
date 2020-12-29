using WriteVTK

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