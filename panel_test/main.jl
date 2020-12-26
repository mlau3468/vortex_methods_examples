include("inf.jl")
include("aeroel.jl")
include("vis.jl")

# solver input
uinf = [30; 0; 0]
Pinf = 0
rhoinf = 1.225
dt = 0.005
te_scaling = 0.3



# Read input file
fname = "./dust_output/wing/geo_input.h5"
fid = h5open(fname, "r")
#display(fid["Components"])


numComp = read(fid["Components"]["NComponents"])
comp_keys = keys(fid["Components"])[1:end-1]
#display(comp_keys)


# Build reference frames
refs = []
ref_keys = []
newRef = refFrame("test", "0", [], false, false, [0;0;0], [0.9961947 0.0 -0.0871557; 0.0 1.0 0.0; 0.0871557 0.0 0.9961947]', [0;0;0], zeros(3,3))
refs = append!(refs, [newRef])
ref_keys = ["test"]

# Build component lists
comps = []
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
    te_dir = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["t_te"])
    n_pan = size(ee,2)
    n_vert = size(rr, 2)
    n_wake_pan = size(ii_te,2)
    n_wake_vert = size(rr_te, 2)

    newComp = component(name, elType, refId, refTag, ee, rr, neigh,  n_pan, n_vert, e_te, i_te, ii_te, rr_te, neigh_te, n_wake_pan, n_wake_vert, te_dir)
    global comps = append!(comps, [newComp])
    global nPan = nPan + n_pan
    global nVert = nVert + n_vert
    global nWakePan = nWakePan + n_wake_pan
    global nWakeVert = nWakeVert + n_wake_vert
end

@time begin
# build global geometry
panels = []
rr_all = zeros(3, nVert)

last_pan = 0
last_vert = 0
for i = 1:numComp
    res = comps[i].refId
    orig = refs[res].origin
    orient = refs[res].orient

    # append panels
    for j = 1:comps[i].n_pan
        ee_comp = comps[i].ee[:,j]
        ee_global = Int.(ee_comp .+ last_vert)
        pts = copy(comps[i].rr[:,ee_comp])
        for k = 1:size(pts,2)
            pts[:,k] = orient*pts[:,k]  .+ orig
        end

        velbody = [0.0; 0.0; 0.0]
        velvort = [0.0; 0.0; 0.0]
        n_ver, n_sides, center, edge_vec, edge_len, edge_uni, tang, normal, area, sinTi, cosTi = getPanProp(pts)
        new_pan = panel(res, i, ee_global, ee_comp, center, n_sides, n_ver, edge_len, edge_vec, edge_uni, normal, area, tang, cosTi, sinTi, velbody, velvort, [0.0])
        push!(panels, new_pan)
        #panels[i] = new_pan
    end
    
    # rotate vertices accordingly and append vertices
    for j = 1:comps[i].n_vert
        rr_all[:, last_vert+1] = orient*comps[i].rr[:,j] + orig
        global last_vert = last_vert + 1
    end

    global last_pan = last_pan + comps[i].n_pan
end
end

@time begin
# set up linear system
A = zeros(nPan, nPan)
B = zeros(nPan, nPan) #Bstatic
RHS = zeros(nPan)
uvort = zeros(3,nPan)
magPan = zeros(nPan)
magWake = zeros(nWakePan)

for i = 1:nPan
    for j = 1:nPan
        if i == j
            dou = -2*pi
        else
            dou = dub(panels[j], panels[i].center, rr_all)
        end
        A[i,j] = -dou

        sou = sourc(panels[j], panels[i].center, dou, rr_all)
        B[i,j] = sou
    end
    RHS[i] = 0.0
end
end

# Build trailing edge
wake_panels = []
rr_wake = zeros(3, nWakeVert)
rr_wake_new = zeros(3, 1)
num_new = 0
last_vert = 0
wake_pan = 0
wake_vert = 0

comp_pan = 0
for i = 1:numComp
    res = comps[i].refId
    ref_vel = refs[res].vel_g
    ref_rot = refs[res].rot_g
    orig = refs[res].origin
    orient = refs[res].orient

    # rotate vertices accordingly and append first row
    for j = 1:comps[i].n_vert_te
        rr_wake[:,last_vert+1] = orient*comps[i].rr_te[:,j] + orig # WARNING: THIS IS IN PLACE
        global last_vert = last_vert + 1
    end

    for j = 1:comps[i].n_pan_te
        ee_comp = comps[i].ii_te[:,j]
        pan_idxcomp = comps[i].e_te[:,j]

        ee_global = ee_comp .+ wake_vert
        pan_idxglobal = pan_idxcomp .+ comp_pan

        pts = copy(comps[i].rr_te[:,ee_comp])
        for k = 1:size(pts,2)
            pts[:,k] = orient*pts[:,k] .+ orig
        end

        # build second row
        w_start_pt = mean(pts, dims=2) # Middle point between these two points
        vel_te = calc_node_vel(w_start_pt, ref_rot, ref_vel)
        te_pan_dir = comps[i].dir_te[:,j]
        dist = orient*te_pan_dir
        new_pts = pts .+ dist.*te_scaling.*norm(uinf.-vel_te).*dt ./ norm(dist)
        global rr_wake_new = [rr_wake_new new_pts]
        ee_global = [ee_global; num_new + nWakeVert + 2]
        ee_global = [ee_global; num_new + nWakeVert + 1]
        global num_new = num_new + 2
        global wake_pan = wake_pan + 1

        # create panel object
        velbody = [0; 0; 0]
        velvort = [0; 0; 0]
        n_ver, n_sides, cpt, edge_vec, edge_len, edge_uni, tang, normal, area, sinTi, cosTi = getPanProp([pts new_pts[:,2] new_pts[:,1]])
        new_pan = wake_panel(res, i, ee_global, ee_comp, pan_idxglobal, pan_idxcomp, te_pan_dir,cpt, n_sides, n_ver, edge_len, edge_vec, edge_uni, normal, area, tang, cosTi, sinTi, velbody, velvort, [0.0])
        push!(wake_panels, new_pan)
    end
    global comp_pan = comp_pan + comps[i].n_pan
    global wake_vert = wake_vert + comps[i].n_vert_te
end
#combine points from first and second row
rr_wake = [rr_wake rr_wake_new[:,2:end]]

panels2vtk(panels, rr_all, "mesh.vtu")
panels2vtk(wake_panels, rr_wake, "mesh_wake.vtu")

# influence of trailing edge
for i = 1:size(panels,1)
    # go through each wake panel
    for j = 1:size(wake_panels,1)
        # effect of wake panel on panel i
        a = dub(wake_panels[j], panels[i].center, rr_wake)
        a = - a # note flipping is down outside comput_pot
        #=Modify influence coefficients of trailing edge panels 
        associated with this trailing edge wake on panel i =#
        A[i, wake_panels[j].panIdx[1]] = A[i, wake_panels[j].panIdx[1]] + a
        A[i, wake_panels[j].panIdx[2]] = A[i, wake_panels[j].panIdx[2]] - a
    end
end

# assemble RHS (add uinf and uvort contribution)
for i = 1:size(panels,1)
    for j = 1:size(panels,1)
        RHS[i] = RHS[i] + B[i,j] .* sum(panels[j].norm.*(-uinf.-panels[j].velVort))
    end
end

# LU decomposition of A
#A = factorize2(A)
# debug: output
writedlm("A.csv", A)
writedlm("B.csv", RHS)
writedlm("B_static.csv", B)


# solve linear system
solution = A\RHS
# debug: output
writedlm("res.csv", solution)

function elemVel(panels, loc, uinf, rr)
    vel_total = zeros(3,0)
    for i = 1:1
    #for i =1:size(panels,1)
       vdub = vel_dub(panels[i], loc, rr)

       vsou = vel_sourc(panels[i], loc, rr)
       nor = panels[i].norm
       mag = panels[i].mag[1]
       uvort = panels[i].velVort
       ub = panels[i].velBody
       vel = vdub.*mag .- vsou .*(sum(nor.*(ub.-uinf.-uvort)))
       if i == 1
       #println(vdub)
       println(vsou)
       end

       vel_total = vel_total.+vel./(4*pi)
    end
    return vel_total
end

# update strengths
for i = 1:size(panels, 1)
    panels[i].mag[1] = solution[i]
end

for i = 1:size(wake_panels,1)
    wake_panels[i].mag[1] = panels[wake_panels[i].panIdx[1]].mag[1] - panels[wake_panels[i].panIdx[2]].mag[1]
end

#for i = 1:size(wake_panels,1)
for i = 1:1
    posp1 = rr_wake[:,wake_panels[i].ee[3]]
    #test = elemVel(panels, posp1, uinf, rr_all)
    
    posp2 = rr_wake[:,wake_panels[i].ee[4]]
    test = elemVel(panels, posp2, uinf, rr_all)

    # left side
    dir1 = rr_wake[:,wake_panels[i].ee[4]]
    
end
