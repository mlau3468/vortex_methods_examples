include("inf.jl")
include("aeroel.jl")
include("vis.jl")
import Base.Threads.@spawn

# solver input
uinf = [30; 0; 0]
Pinf = 0
rhoinf = 1.225
dt = 0.005
te_scaling = 0.3
t_end = 100*dt
debug_dir = "./debug/"
vis_dir = "./vis/"

t = 0.0
step_num = Int(1)



# Read input file
fname = "./dust_output/wing/geo_input.h5"
fid = h5open(fname, "r")

#display(fid["Components"])
#display(read(fid["Components"]["Comp001"]["Trailing_Edge"]["o_te"]))


numComp = read(fid["Components"]["NComponents"])
comp_keys = keys(fid["Components"])[1:end-1]
#display(comp_keys)


# Build reference frames
refs = []
ref_keys = []
newRef = refFrame("test", "0", [], false, false, [0;0;0], [0.9961947 0.0 -0.0871557; 0.0 1.0 0.0; 0.0871557 0.0 0.9961947]', [0;0;0], zeros(3,3))
# orient is defined as v_global = orient * v_local
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
    neigh_te_o = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["o_te"])
    te_dir = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["t_te"])
    n_pan = size(ee,2)
    n_vert = size(rr, 2)
    n_wake_pan = size(ii_te,2)
    n_wake_vert = size(rr_te, 2)

    newComp = component(name, elType, refId, refTag, ee, rr, neigh,  n_pan, n_vert, e_te, i_te, ii_te, rr_te, neigh_te, neigh_te_o, n_wake_pan, n_wake_vert, te_dir)
    global comps = append!(comps, [newComp])
    global nPan = nPan + n_pan
    global nVert = nVert + n_vert
    global nWakePan = nWakePan + n_wake_pan
    global nWakeVert = nWakeVert + n_wake_vert
end

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

# Build trailing edge and prep wake
wake_panels = []
wake_particles = []
end_vorts = []
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
        if pan_neighglobal[1] > 0
            pan_neighglobal[1] = pan_neighglobal[1] .+ wake_pan
        end
        if pan_neighglobal[2] > 0
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
        global rr_wake_new = [rr_wake_new new_pts]
        ee_global = [ee_global; num_new + nWakeVert + 2]
        ee_global = [ee_global; num_new + nWakeVert + 1]
        global num_new = num_new + 2

        # create panel object
        velbody = [0; 0; 0]
        velvort = [0; 0; 0]
        n_ver, n_sides, cpt, edge_vec, edge_len, edge_uni, tang, normal, area, sinTi, cosTi = getPanProp([pts new_pts[:,2] new_pts[:,1]])
        new_pan = wake_panel(res, i, ee_global, ee_comp, pan_idxglobal, pan_idxcomp, te_pan_dir,cpt, n_sides, n_ver, edge_len, edge_vec, edge_uni, normal, area, tang, cosTi, sinTi, velbody, velvort, [0.0], pan_neighglobal, pan_neighcomp, neigh_te_orient)
        push!(wake_panels, new_pan)
    end
    global comp_pan = comp_pan + comps[i].n_pan
    global wake_vert = wake_vert + comps[i].n_vert_te
    global wake_pan = wake_pan + comps[i].n_pan_te
end
#combine points from first and second row
rr_wake = [rr_wake rr_wake_new[:,2:end]]

lastpanidou = zeros(nWakePan) #Last vortex intensity from removed panels
endpanidou = zeros(nWakePan) #Last vortex intensity from removed panels

# visualize
panels2vtk(panels, rr_all, "mesh_$step_num.vtu", vis_dir)
panels2vtk(wake_panels, rr_wake, "wake_pan_$step_num.vtu", vis_dir)
particles2vtk(wake_particles, "wake_particles_$step_num.vtu", vis_dir)
step_num = step_num + 1














# Timestepping
while t < t_end
    # update uvort, influence of vortex particles and end vortex line

for i = 1:size(panels, 1)
        vel = uvortParticles(wake_particles, panels[i].center)
        vel2 = uvortVortLines(end_vorts, panels[i].center)
        vel = vel .+ vel2
        panels[i].velVort[:] = vel[:]
        #print(i)
        #print(' ')
        #println(panels[i].velVort[:])
    end
#=
    for i = 1:size(wake_panels,1)
        vel = uvortParticles(wake_particles, panels[i].center)
        vel2 = uvortVortLines(end_vorts, panels[i].center)
        vel = vel .+ vel2
        wake_panels[i].velVort[:] = vel[:]
    end
=#  

# influence of panels
# set up linear system
A = zeros(nPan, nPan)
B = zeros(nPan, nPan) #Bstatic
RHS = zeros(nPan)

Threads.@threads for i = 1:nPan
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

# influence of trailing edge
for i = 1:size(panels,1)
    # go through each wake panel
    for j = 1:size(wake_panels,1)
        # effect of wake panel on panel i
        local a = dub(wake_panels[j], panels[i].center, rr_wake)
        local a = - a # note flipping is down outside comput_pot
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


# debug: output
writedlm("B_static.csv", B)

# solve linear system
solution = A\RHS
# debug: output
debugMatrices(A, RHS, solution, step_num-1, debug_dir)

# update strengths
for i = 1:size(panels, 1)
    panels[i].mag[1] = solution[i]
end
for i = 1:size(wake_panels,1)
    wake_panels[i].mag[1] = panels[wake_panels[i].panIdx[1]].mag[1] - panels[wake_panels[i].panIdx[2]].mag[1]
    #println(wake_panels[i].mag[1])
end

# calculate velocities at existing particles
for i = 1:size(wake_particles,1)
    vel = elemVel(panels, wake_panels, wake_particles, end_vorts, rr_all, rr_wake, wake_particles[i].center) .+ uinf
    wake_particles[i].vel[:] = vel
end

# update particle positions
for i = 1:size(wake_particles,1)
    wake_particles[i].center[:] = wake_particles[i].center .+ wake_particles[i].vel .* dt
    #println(wake_particles[i].center[:])
end

# shed new particles
pts1 = zeros(3,nWakePan) # temporary vector to hold end vortex points
pts2 = zeros(3,nWakePan)
for i = 1:size(wake_panels,1)
    posp1 = rr_wake[:,wake_panels[i].ee[3]]
    v1 = elemVel(panels, wake_panels, wake_particles, end_vorts, rr_all, rr_wake, posp1) .+ uinf
    posp2 = rr_wake[:,wake_panels[i].ee[4]]
    v2 = elemVel(panels, wake_panels, wake_particles, end_vorts, rr_all, rr_wake, posp2) .+ uinf
    pts1[:,i] = posp1 .+ v1.*dt
    pts2[:,i] = posp2 .+ v2.*dt

    #println(posp1)
    #println(pts1[:,i])
    #println(posp2)
    #println(pts1[:,i])
    #println("------")

end

for i = 1:size(wake_panels,1)
    partvec = zeros(3)
    # left side
    #println(wake_panels[i].neigh_te)
    dir = rr_wake[:,wake_panels[i].ee[4]] .- pts2[:,i]
    if wake_panels[i].neigh_te[1] > 0
        ave = wake_panels[i].mag[1] - wake_panels[i].neigh_orient[1] * wake_panels[wake_panels[i].neigh_te[1]].mag[1]
        ave = ave/2
    else
        ave = wake_panels[i].mag[1]
    end
    partvec = partvec .+ dir.*ave

    # right side
    dir = -rr_wake[:,wake_panels[i].ee[3]] .+ pts1[:,i]
    if wake_panels[i].neigh_te[2] > 0
        ave = wake_panels[i].mag[1] - wake_panels[i].neigh_orient[2] * wake_panels[wake_panels[i].neigh_te[2]].mag[1]
        ave = ave/2
    else
        ave = wake_panels[i].mag[1]
    end
    partvec = partvec .+ dir.*ave

    # end side
    dir =  pts2[:,i] .- pts1[:,i]
    ave = wake_panels[i].mag[1]-lastpanidou[i]
    lastpanidou[i] = wake_panels[i].mag[1]
    partvec = partvec .+ dir.*ave
    
    # calculate the center
    posp = (pts1[:,i] .+ pts2[:,i] .+ rr_wake[:,wake_panels[i].ee[4]].+ rr_wake[:,wake_panels[i].ee[3]])./4

    # add wake particle
    cent = posp
    mag = norm(partvec)
    if wake_panels[i].mag[1] > 1e-13
        dir = partvec./mag
    else
        dir = partvec
    end
    vel = [0;0;0] # panel doesn't move
    new_part = wake_part(dir, [mag], cent, vel)
    #println(cent)
    push!(wake_particles, new_part)
    endpanidou[i] = wake_panels[i].mag[1]

    # attach end vortex
    ee = [wake_panels[i].ee[4]; wake_panels[i].ee[3]]
    p1 = rr_wake[:,wake_panels[i].ee[4]]
    p2 = rr_wake[:,wake_panels[i].ee[3]]
    mag = lastpanidou[i]
    mag = [mag]
    ver_vel = [0;0;0]
    center, n_ver, n_sides, edge_vec, edge_len, edge_uni = getLineProp([p1 p2])
    newvortline = vortline(ee,  [p1 p2], center, edge_vec, edge_uni, edge_len, ver_vel, mag)

    if size(end_vorts,1) < size(wake_panels,1)
        push!(end_vorts, newvortline)
    else
        end_vorts[i] = newvortline
    end

end

# visualize
panels2vtk(panels, rr_all, "mesh_$step_num.vtu", vis_dir)
panels2vtk(wake_panels, rr_wake, "wake_pan_$step_num.vtu", vis_dir)
particles2vtk(wake_particles, "wake_particles_$step_num.vtu", vis_dir)

# attach end vortex

global t = t + dt
global step_num = step_num + Int(1)

println("Time: $t")
end # timestepping