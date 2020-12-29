struct panel
    refId :: Int
    compId :: Int
    ee  :: Array{Int32,1}# index of points forming vertices in global points list
    eeComp :: Array{Int32,1}# index of points forming vertices in local points list
    center:: Array{Float64,1}
    nSide :: Int
    nVert :: Int
    edgeLen ::Array{Float64,1}
    edgeVec ::Array{Float64,2}
    edgeUni ::Array{Float64,2}
    norm ::Array{Float64,1} # element unit normal vector
    area :: Float64
    tang ::Array{Float64,2} # teangent unit vectors as in PANAIR
    cosTi ::Array{Float64,1}
    sinTi ::Array{Float64,1}
    velBody ::Array{Float64,1} # body velocity at the center, global coordinates
    velVort ::Array{Float64,1} # vorticity induced velocity at the center
    mag :: Array{Float64,1}
end

struct wake_panel
    refId :: Int
    compId :: Int
    ee  :: Array{Int32,1}# index of points forming vertices in global points list
    eeComp :: Array{Int32,1}# index of points forming vertices in local points list
    panIdx :: Array{Int32,1} # index of panel trailing edge is attached to (global)
    pandIdxComp :: Array{Int32,1} # index of panel trailing edge is attached to (component)
    dirTE :: Array{Float64,1} #Unit vector at TE nodes
    center:: Array{Float64,1}
    nSide :: Int
    nVert :: Int
    edgeLen ::Array{Float64,1}
    edgeVec ::Array{Float64,2}
    edgeUni ::Array{Float64,2}
    norm ::Array{Float64,1} # element unit normal vector
    area :: Float64
    tang ::Array{Float64,2} # teangent unit vectors as in PANAIR
    cosTi ::Array{Float64,1}
    sinTi ::Array{Float64,1}
    velBody ::Array{Float64,1} # body velocity at the center
    velVort ::Array{Float64,1} # vorticity induced velocity at the center
    mag :: Array{Float64,1}
    neigh_te :: Array{Int,1} # neighboring wake element indices global
    neigh_te_comp :: Array{Int,1} # neighboring wake element indices in component
    neigh_orient :: Array{Float64,1} # relative orientation of each neighboring element
end

struct wake_part
    dir :: Array{Float64, 1}
    mag :: Array{Float64,1}
    center :: Array{Float64,1}
    vel :: Array{Float64,1}
end

function getPanProp(pts)
    n_ver = 4
    n_sides = 4

    # settings
    prev_qua = [4 1 2 3]
    next_qua = [2 3 4 1]
    prev_tri = [3 1 2]
    next_tri = [2 3 1]

    # get edge vectors
    edge_vec = zeros(3,4)

    if n_sides == 3
    elseif n_sides == 4
        for i = 1:n_sides
            edge_vec[:,i] = pts[:, next_qua[i]] .- pts[:, i]
        end
    end

    # get edge lens
    edge_len = zeros(4)
    # unit vectors
    edge_uni = zeros(3,4)
    for i = 1:n_sides
        edge_len[i] = norm(edge_vec[:,i])
        edge_uni[:,i] = edge_vec[:,i] ./ norm(edge_vec[:,i])
    end

    # normal 
    v3 = cross(edge_vec[:,1], edge_vec[:,2])
    normal = v3./norm(v3)

    # central point
    cpt = mean(pts, dims=2)[:,1]

    # local tangent unit vector as in PANAIR
    tanl = 0.5 .* (pts[:,n_sides] .+ pts[:,1]) .- cpt
    tang = zeros(3,2)
    tang[:,1] = tanl ./ norm(tanl)
    tang[:,2] = cross(normal, tang[:,1])



    # area
    area = 0.5*norm(v3)

    sinTi = zeros(n_sides)
    cosTi = zeros(n_sides)

    for i = 1:n_sides
        cosTi[i] = sum(edge_uni[:,i] .* tang[:,1])
        sinTi[i] = sum(edge_uni[:,i] .* tang[:,2])
    end

    return(n_ver, n_sides, cpt, edge_vec, edge_len, edge_uni, tang, normal, area, sinTi, cosTi)
end

struct refFrame
    name :: String
    parent :: String
    children # list of child references
    move2Parent :: Bool # Is this reference frame moving with respect to the parent frame
    move2Global :: Bool # Is this reference frame moving with respect to the global frame
    origin # 3x1 vector origin in the parent frame
    orient # 3x3 orientation matrix with respect to parent

    vel_g # 3x1 vector velocity with respect to global frame 
    rot_g #3x3 rotation rate matrix with respect to global frame
end

struct component
    name::String
    elType :: String
    refId :: Int
    refName :: String
    ee  # panels
    rr  # vertices
    neigh # neighboring
    n_pan :: Int # Number of panels
    n_vert :: Int # Number of vertices

    e_te #panel index associated with each TE panel
    i_te #component panel id of the nodes on the TE
    ii_te #TE id of the nodes of the TE elements (index in rr_te)
    rr_te #Coordinates of the nodes of the TE
    neigh_te #TE id of neighboring TE elements
    neigh_te_o # relative orientation of neighbors
    n_pan_te :: Int
    n_vert_te :: Int
    dir_te #Unit vector at TE nodes

end 