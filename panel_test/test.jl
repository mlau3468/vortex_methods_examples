include("panel.jl")
using WriteVTK

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

# solver input
uinf = [30; 0; 0]
Pinf = 0
rhoinf = 1.225
dt = 0.005

# global points

ee_all = zeros(4, nPan)
rr_all = zeros(3, nVert)

last_pan = 0
last_vert = 0
for i = 1:numComp
    res = comps[i].refId
    orig = refs[res].origin
    orient = refs[res].orient

    # append panels
    ee_all[ :, last_pan+1:last_pan+comps[i].n_pan] = comps[i].ee .+ last_vert
    
    # rotate vertices accordingly and append vertices
    for j = 1:comps[i].n_vert
        rr_all[:, last_vert+1] = orient*comps[i].rr[:,j] + orig
        global last_vert = last_vert + 1
    end

    global last_pan = last_pan + comps[i].n_pan
end

# recast as integers
ee_all = Int.(ee_all)

# center points of each panel
cent = zeros(3, size(ee_all,2))
for j = 1:size(ee_all,2)
    cent[:, j] = mean(rr_all[:, ee_all[:, j]], dims=2)
end

# calculate normals
normals = zeros(3, size(ee_all,2))
for j = 1:size(ee_all,2)
    v1 = rr_all[:, ee_all[2, j]] .- rr_all[:, ee_all[1, j]]
    v2 = rr_all[:, ee_all[3, j]] .- rr_all[:, ee_all[2, j]]
    v3 = cross(v1, v2)
    normals[:, j] = v3./norm(v3)
end

# set up linear system
A = zeros(nPan, nPan)
B = zeros(nPan, nPan)
for i = 1:nPan
    for j = 1:nPan
        if i == j
            dou = -2*pi
        else
            dou = dub(rr_all[:,ee_all[:,j]], cent[:,i], normals[:,j])
        end
        A[i,j] = -dou

        sou = sourc(rr_all[:,ee_all[:,j]], cent[:,i], dou, normals[:,j])
        B[i,j] = sou
    end
end

# trailing edge first row
last_pan = 0
last_vert = 0
comp_pan  = 0
ee_wake = zeros(4, nWakePan)
rr_wake = zeros(3, nWakeVert)
te_pan = zeros(2, nWakePan)
te_scaling = 0.3

for i = 1:numComp
    res = comps[i].refId
    ref_vel = refs[res].vel_g
    ref_rot = refs[res].rot_g
    orig = refs[res].origin
    orient = refs[res].orient

    ee_wake[ 1:2, last_pan+1:last_pan+comps[i].n_pan_te] = comps[i].ii_te .+ last_vert
    # rotate vertices accordingly
    for j = 1:comps[i].n_vert_te
        rr_wake[:,last_vert+1] = orient*comps[i].rr_te[:,j] + orig # WARNING: THIS IS IN PLACE
        global last_vert = last_vert + 1
    end

    # figure out which panels are associated with the trailing edge
    for j = 1:comps[i].n_pan_te
        te_pan[:,last_pan.+j] = comp_pan .+ comps[i].e_te[:,j]
    end

    global last_pan = last_pan + comps[i].n_pan_te
    global comp_pan = comp_pan + comps[i].n_pan
end


# trailing edge second row
last_pan = 0
for i = 1:numComp
    res = comps[i].refId
    ref_vel = refs[res].vel_g
    ref_rot = refs[res].rot_g
    orig = refs[res].origin
    orient = refs[res].orient
    for j = 1:comps[i].n_pan_te
        pts = comps[i].rr_te[:,comps[i].ii_te[:,j]]
        pts = [orient*pts[:,1] + orig  orient*pts[:,2] + orig] # coordinate frame transformation, cancel due to not copying
        w_start_pt = mean(pts, dims=2) # Middle point between these two points
        vel_te = calc_node_vel(w_start_pt, ref_rot, ref_vel)
        te_pan_dir = comps[i].dir_te[:,j]
        dist = orient*te_pan_dir
        new_pts = pts .+ dist.*te_scaling.*norm(uinf.-vel_te).*dt ./ norm(dist)
        global rr_wake = [rr_wake new_pts] # append
        ee_wake[3:4, last_pan+1] = [last_vert+2;last_vert+1]
        global last_pan = last_pan + 1
        global last_vert = last_vert + 2
    end  
end


ee_wake = Int.(ee_wake)
te_pan = Int.(te_pan)


# influence of trailing edge
for i = 1:size(ee_all,2)
    # go through each wake panel
    for j = 1:size(ee_wake,2)
        # effect of wake panel on panel i
        a = dub(rr_wake[:,ee_wake[:,j]], cent[:,i], normals[:,j])
        a = - a # note flipping is down outside comput_pot
        #=Modify influence coefficients of trailing edge panels 
        associated with this trailing edge wake on panel i =#
        A[i,te_pan[1,j]] = A[i,te_pan[1,j]] + a
        A[i,te_pan[2,j]] = A[i,te_pan[2,j]] - a
    end
    
    
end

# LU decomposition of A
A = factorize2(A)
writedlm("A.csv", A)
writedlm("B.csv", B)

writedlm("ee.csv", ee_wake', ',')
writedlm("rr.csv", rr_wake', ',')
eerr2vtk(ee_all, rr_all, "mesh.vtu")

#=
C = readdlm("test.txt", ',')
D = zeros(nPan, nPan)
for i=1:size(C,1)
    D[Int(C[i,1]), Int(C[i,2])] = C[i,3]
end
=#
#display(fid["Components"])
#=
display(fid["Components"]["Comp001"]["Trailing_Edge"]) 
display(read(fid["Components"]["Comp001"]["Trailing_Edge"]["e_te"])) #Global id of the elements at the TE
display(read(fid["Components"]["Comp001"]["Trailing_Edge"]["i_te"])) #Global id of the nodes on the TE
display(read(fid["Components"]["Comp001"]["Trailing_Edge"]["ii_te"])) #TE id of the nodes of the TE elements
display(read(fid["Components"]["Comp001"]["Trailing_Edge"]["neigh_te"])) #TE id of neighboring TE elements
display(read(fid["Components"]["Comp001"]["Trailing_Edge"]["o_te"])) #Relative orientation of the neighboring TE elements
display(read(fid["Components"]["Comp001"]["Trailing_Edge"]["rr_te"])) #Coordinates of the nodes of the TE
display(read(fid["Components"]["Comp001"]["Trailing_Edge"]["scale_te"])) #Individual scaling for each component
display(read(fid["Components"]["Comp001"]["Trailing_Edge"]["t_te"])) #Unit vector at TE nodes
=#


