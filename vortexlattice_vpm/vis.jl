using WriteVTK

function particles2vtk(particles_list, fname)

    npoints = size(particles_list, 1)
    if npoints > 0
        x = zeros(size(particles_list, 1))
        y = zeros(size(particles_list, 1))
        z = zeros(size(particles_list, 1))
        mag = zeros(size(particles_list, 1))
        for i = 1:size(particles_list, 1)
            x[i] = particles_list[i].cpt[1]
            y[i] = particles_list[i].cpt[2]
            z[i] = particles_list[i].cpt[3]
            mag[i] = particles_list[i].gam[1]
        end
    else
        npoints = 1
        x = [NaN]
        y = [NaN]
        z = [NaN]
        mag = [NaN]
    end

    cells = [MeshCell(VTKCellTypes.VTK_VERTEX, (i, )) for i = 1:npoints]

    vtk_grid(fname, x, y, z, cells) do vtk
        vtk["mag", VTKPointData()] = mag
    end
end

function panels2vtk(panel_list, fname)
    pts = zeros(3,length(panel_list)*4)
    for i = 1:length(panel_list)
        for j=1:4
            pts[:,(i-1)*4+j] = panel_list[i].pts[:,j]
        end
    end

    celltype = VTKCellTypes.VTK_QUAD
    cells = MeshCell[]
    mag = zeros(size(panel_list,1))
    #velmag = zeros(size(panel_list,1))
    pres = zeros(size(panel_list,1))
    inds = [1;2;3;4]
    for i =1:size(panel_list,1)
        c = MeshCell(celltype, inds)
        push!(cells, c)
        mag[i] = panel_list[i].gam[1]
        pres[i] = panel_list[i].dp[1]
        inds = inds .+ 4
    end
    outfile = vtk_grid(fname, pts, cells, compress=2) do vtk
        vtk["mag"] = mag
        #vtk["velmag"] = velmag
        vtk["pres"] = pres
    end
end

function wakepanels2vtk(panel_list, fname)
    pts = zeros(3,length(panel_list)*4)
    for i = 1:length(panel_list)
        for j=1:4
            pts[:,(i-1)*4+j] = panel_list[i].pts[:,j]
        end
    end

    celltype = VTKCellTypes.VTK_QUAD
    cells = MeshCell[]
    inds = [1;2;3;4]
    for i =1:size(panel_list,1)
        c = MeshCell(celltype, inds)
        push!(cells, c)
        inds = inds .+ 4
    end
    outfile = vtk_grid(fname, pts, cells, compress=2) do vtk
    end
end
