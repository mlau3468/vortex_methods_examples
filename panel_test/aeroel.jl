struct panel
    refId :: Int
    compId :: Int
    ee  :: Array{Int32,1}# index of points forming vertices in global points list
    eeComp :: Array{Int32,1}# index of points forming vertices in local points list
    center:: Array{Float64,1}
    neigh :: Array{Int32,1} # index of neighboring elements in global
    neighComp :: Array{Int32,1} # index of neighboring elements in component
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
    mag :: Array{Float64,1} # magnitude
    center :: Array{Float64,1} # position in global frame
    vel :: Array{Float64,1} # velocity in global frame
end

struct vortline
    ee :: Array{Int32,1} # index in rr wake
    rr :: Array{Float64,2} #points
    center ::Array{Float64,1}
    edgeVec ::Array{Float64,1}
    edgeUni ::Array{Float64,1}
    edgeLen :: Float64
    vel :: Array{Float64,1}
    mag :: Array{Float64,1} # magnitude
end 

function getLineProp(pts)
    center =  mean(pts, dims=2)[:,1]
    n_ver = 2
    n_sides = 1

    edge_vec = pts[:,2] .- pts[:,1]
    edge_len = norm(edge_vec)
    edge_uni = edge_vec ./ edge_len
    return(center, n_ver, n_sides, edge_vec, edge_len, edge_uni)
end

function getPanProp(pts)

    n_ver = size(pts,2)
    n_sides = size(pts, 2)

    # settings
    prev_qua = [4 1 2 3]
    next_qua = [2 3 4 1]
    prev_tri = [3 1 2]
    next_tri = [2 3 1]

    # get edge vectors
    edge_vec = zeros(3,n_sides)

    if n_sides == 3
        for i = 1:n_sides
            edge_vec[:,i] = pts[:, next_tri[i]] .- pts[:, i]
        end
    elseif n_sides == 4
        for i = 1:n_sides
            edge_vec[:,i] = pts[:, next_qua[i]] .- pts[:, i]
        end
    end

    # get edge lens
    edge_len = zeros(n_sides)
    # unit vectors
    edge_uni = zeros(3,n_sides)
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

function get_geo(fname, uinf, refs, ref_keys)
te_scaling = 0.3
# Read input file
#fname = "./dust_output/wing/geo_input.h5"
fid = h5open(fname, "r")

#display(fid["Components"])
#display(read(fid["Components"]["Comp001"]["Trailing_Edge"]["o_te"]))


numComp = read(fid["Components"]["NComponents"])
comp_keys = keys(fid["Components"])[1:end-1]
#display(comp_keys)


# Build component lists
comps = component[]
nPan = 0
nVert = 0
nWakePan = 0
nWakeVert = 0

for i = 1:numComp
    name = read(fid["Components"][comp_keys[i]]["CompName"])
    ee = read(fid["Components"][comp_keys[i]]["Geometry"]["ee"])
    neigh = read(fid["Components"][comp_keys[i]]["Geometry"]["neigh"])
    rr = read(fid["Components"][comp_keys[i]]["Geometry"]["rr"])
    refTag = read(fid["Components"][comp_keys[i]]["RefTag"])
    refId = findall(x -> x.==refTag, ref_keys)[1]
    elType = read(fid["Components"][comp_keys[i]]["ElType"])

    e_te = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["e_te"])
    i_te = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["i_te"])
    ii_te = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["ii_te"])
    rr_te = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["rr_te"])
    neigh_te = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["neigh_te"])
    neigh_te_o = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["o_te"])
    te_dir = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["t_te"])
    n_pan = size(ee,2)
    n_vert = size(rr, 2)
    n_wake_pan = size(ii_te,2)
    n_wake_vert = size(rr_te, 2)

    newComp = component(name, elType, refId, refTag, ee, rr, neigh,  n_pan, n_vert, e_te, i_te, ii_te, rr_te, neigh_te, neigh_te_o, n_wake_pan, n_wake_vert, te_dir)
    comps = append!(comps, [newComp])
    nPan = nPan + n_pan
    nVert = nVert + n_vert
    nWakePan = nWakePan + n_wake_pan
    nWakeVert = nWakeVert + n_wake_vert
end

# build global geometry
panels = panel[]
rr_all = zeros(3, nVert)

last_pan = 0
last_vert = 0
for i = 1:numComp
    res = comps[i].refId
    orig = refs[res].origin
    orient = refs[res].orient

    # append panels
    for j = 1:comps[i].n_pan
        # local component indexes
        ee_comp = comps[i].ee[:,j]
        neigh_comp = comps[i].neigh[:,j]
    
        # if triangle, last vertex is filled with a 0
        if ee_comp[end] == 0
            ee_comp = ee_comp[1:end-1]
            neigh_comp = neigh_comp[1:end-1]
        end

        # global indexes
        ee_global = Int.(ee_comp .+ last_vert)

        neigh_global =  copy(neigh_comp)
        for k = 1:size(neigh_comp,1) # ignore 0s, as they mean no neighbor on this edge
            if neigh_global[k] != 0
            neigh_global[k] = Int.(neigh_global[k] .+ last_pan)
            end
        end
        
        pts = copy(comps[i].rr[:,ee_comp])
        for k = 1:size(pts,2)
            pts[:,k] = orient*pts[:,k]  .+ orig
        end

        velbody = [0.0; 0.0; 0.0]
        velvort = [0.0; 0.0; 0.0]
        n_ver, n_sides, center, edge_vec, edge_len, edge_uni, tang, normal, area, sinTi, cosTi = getPanProp(pts)
        new_pan = panel(res, i, ee_global, ee_comp, center, neigh_global, neigh_comp, n_sides, n_ver, edge_len, edge_vec, edge_uni, normal, area, tang, cosTi, sinTi, velbody, velvort, [0.0])
        push!(panels, new_pan)
        #panels[i] = new_pan
    end
    
    # rotate vertices accordingly and append vertices
    for j = 1:comps[i].n_vert
        rr_all[:, last_vert+1] = orient*comps[i].rr[:,j] + orig
        last_vert = last_vert + 1
    end

    last_pan = last_pan + comps[i].n_pan
end

# Build trailing edge and prep wake
wake_panels = wake_panel[]
wake_particles = wake_part[]
end_vorts = vortline[]
rr_wake = zeros(3, nWakeVert)
rr_wake_new = zeros(3, 1)
num_new = 0 # counter: number of additional wake vertices added due to second row
wake_vert = 0 # counter: number of wake vertices 
wake_pan = 0 # counter: number of wake panels added
comp_pan = 0 # counter: number of component panels traversed
for i = 1:numComp
    res = comps[i].refId
    ref_vel = refs[res].vel_g
    ref_rot = refs[res].rot_g
    orig = refs[res].origin
    orient = refs[res].orient

    # rotate vertices accordingly and append first row
    for j = 1:comps[i].n_vert_te
        rr_wake[:,wake_vert+j] = orient*comps[i].rr_te[:,j] + orig
    end

    for j = 1:comps[i].n_pan_te
        # wake panels
        ee_comp = comps[i].ii_te[:,j]
        ee_global = ee_comp .+ wake_vert

        # panel connnected to this wake panel
        pan_idxcomp = comps[i].e_te[:,j]
        pan_idxglobal = pan_idxcomp .+ comp_pan

        # panel neighbors
        pan_neighcomp = comps[i].neigh_te[:,j]
        pan_neighglobal = copy(pan_neighcomp)
        if pan_neighglobal[1] > 0 # ignore 0s, as they mean no neighbor on this edge
            pan_neighglobal[1] = pan_neighglobal[1] .+ wake_pan
        end
        if pan_neighglobal[2] > 0 # ignore 0s, as they mean no neighbor on this edge
            pan_neighglobal[2] = pan_neighglobal[2] .+ wake_pan
        end
        pan_neighcomp = Int.(pan_neighcomp)
        pan_neighglobal = Int.(pan_neighglobal)

        # neighbor orientations
        neigh_te_orient = comps[i].neigh_te_o[:,j]

        pts = copy(comps[i].rr_te[:,ee_comp])
        for k = 1:size(pts,2)
            pts[:,k] = orient*pts[:,k] .+ orig
        end

        # build second row
        w_start_pt = mean(pts, dims=2) # Middle point between these two points
        vel_te = calc_node_vel(w_start_pt, ref_rot, ref_vel)
        te_pan_dir = orient*comps[i].dir_te[:,j]
        dist = te_pan_dir
        new_pts = pts .+ dist.*te_scaling.*norm(uinf.-vel_te).*dt ./ norm(dist)
        #println(pts')
        #println(new_pts')
        #println("--------")
        rr_wake_new = [rr_wake_new new_pts]
        ee_global = [ee_global; num_new + nWakeVert + 2]
        ee_global = [ee_global; num_new + nWakeVert + 1]
        num_new = num_new + 2

        # create panel object
        velbody = [0; 0; 0]
        velvort = [0; 0; 0]
        n_ver, n_sides, cpt, edge_vec, edge_len, edge_uni, tang, normal, area, sinTi, cosTi = getPanProp([pts new_pts[:,2] new_pts[:,1]])
        new_pan = wake_panel(res, i, ee_global, ee_comp, pan_idxglobal, pan_idxcomp, te_pan_dir,cpt, n_sides, n_ver, edge_len, edge_vec, edge_uni, normal, area, tang, cosTi, sinTi, velbody, velvort, [0.0], pan_neighglobal, pan_neighcomp, neigh_te_orient)
        push!(wake_panels, new_pan)
    end
    comp_pan = comp_pan + comps[i].n_pan
    wake_vert = wake_vert + comps[i].n_vert_te
    wake_pan = wake_pan + comps[i].n_pan_te
end
#combine points from first and second row
rr_wake = [rr_wake rr_wake_new[:,2:end]]
return panels, rr_all, wake_panels, rr_wake, wake_particles, end_vorts
end

function update_particles(panels, wake_panels, wake_particles, end_vorts, rr_all, rr_wake, uinf, dt)
    # calculate velocities at existing particles
    @Threads.threads for i in 1:size(wake_particles,1)
        vel =  elemVel(panels, wake_panels, wake_particles, end_vorts, rr_all, rr_wake, wake_particles[i].center) .+ uinf
        wake_particles[i].vel[:] = vel
    end
    # update particle positions
    @Threads.threads for i in 1:size(wake_particles,1)
        wake_particles[i].center[:] = wake_particles[i].center .+wake_particles[i].vel.*dt
    end
    return wake_particles
end

function update_panel_uvorts(panels, wake_particles, end_vorts)
    @inbounds @Threads.threads for i = 1:size(panels, 1)
        vel = uvortParticles(wake_particles, panels[i].center)
        vel2 = uvortVortLines(end_vorts, panels[i].center)
        vel = vel .+ vel2
        panels[i].velVort[:] = vel[:]
        #print(i)
        #print(' ')
        #println(panels[i].velVort[:])
    end
    return panels
end