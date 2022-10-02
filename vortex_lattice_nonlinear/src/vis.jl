
function vis_mesh_temp(pan_con::Array{<:Integer}, pan_vert::Array{<:Real}, gam, vmags, fname::String)
    npan = size(pan_con, 2)
    cells = MeshCell[]
    @views for i =1:npan
        inds = pan_con[:,i]
        if inds[end] == 0 # triangle
            c = MeshCell(VTKCellTypes.VTK_TRIANGLE, inds[1:end-1])
        else # quad
            c = MeshCell(VTKCellTypes.VTK_QUAD, inds)
        end
        push!(cells, c)
    end
    vtk_grid(fname, pan_vert, cells, compress=2) do vtk
        vtk["gam"] = gam
        vtk["vmag"] = vmags
    end
end