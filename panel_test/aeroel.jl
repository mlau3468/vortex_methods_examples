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
    velBody ::Array{Float64,1} # body velocity at the center
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
end

function getPanProp(pts)
    n_ver = 4
    n_sides = 4

    # get edge vectors
    edge_vec = zeros(3,4)
    edge_vec[:,1] = pts[:,2] .- pts[:,1]
    edge_vec[:,2] = pts[:,3] .- pts[:,2]
    edge_vec[:,3] = pts[:,4] .- pts[:,3]
    edge_vec[:,4] = pts[:,1] .- pts[:,4]
    # get edge lens
    edge_len = zeros(4)
    edge_len[1] = norm(pts[:,2] .- pts[:,1])
    edge_len[2] = norm(pts[:,3] .- pts[:,2])
    edge_len[3] = norm(pts[:,4] .- pts[:,3])
    edge_len[4] = norm(pts[:,1] .- pts[:,4])

    # unit vectors
    edge_uni = zeros(3,4)
    edge_uni[:,1] = edge_vec[:,1] ./ norm(edge_vec[:,1])
    edge_uni[:,2] = edge_vec[:,2] ./ norm(edge_vec[:,2])
    edge_uni[:,3] = edge_vec[:,3] ./ norm(edge_vec[:,3])
    edge_uni[:,4] = edge_vec[:,4] ./ norm(edge_vec[:,4])

    # central point
    cpt = mean(pts, dims=2)[:,1]
    # normal 
    v3 = cross(edge_vec[:,1], edge_vec[:,2])
    normal = v3./norm(v3)

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
        sinTi[i] = sum(edge_uni[:,i] .* tang[:,1])
        cosTi[i] = sum(edge_uni[:,i] .* tang[:,2])
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
    n_pan_te :: Int
    n_vert_te :: Int
    dir_te #Unit vector at TE nodes

end 