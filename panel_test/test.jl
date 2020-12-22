using HDF5
using Statistics
using LinearAlgebra


fname = "./dust_output/wing/geo_input.h5"
fid = h5open(fname, "r")
#display(fid["Components"])


numComp = read(fid["Components"]["NComponents"])
comp_keys = keys(fid["Components"])[1:end-1]
#display(comp_keys)

struct refFrame
    name :: String
    parent :: String
    children # list of child references
    move2Parent :: Bool # Is this reference frame moving with respect to the parent frame
    move2Global :: Bool # Is this reference frame moving with respect to the global frame
    origin # 3x1 vector rigin in the parent frame
    orient # 3x3 orientation matrix with respect to parent
end

struct component
    name::String
    elType :: String
    refName :: String
    ee  # points
    rr  # vertices
    neigh # neighboring
    ee_te
    rr_te
    neigh_te
    n_pan :: Int
    n_vert :: Int
end 

# Build reference frames
refs = []
ref_keys = []
newRef = refFrame("test", "0", [], false, false, [0;0;0], [0.9961947 0.0 -0.0871557; 0.0 1.0 0.0; 0.0871557 0.0 0.9961947]')
refs = append!(refs, [newRef])
ref_keys = ["test"]

comps = []
nPan = 0
nVert = 0

for i = 1:numComp
    name = read(fid["Components"][comp_keys[i]]["CompName"])
    ee = read(fid["Components"][comp_keys[i]]["Geometry"]["ee"])
    neigh = read(fid["Components"][comp_keys[i]]["Geometry"]["neigh"])
    rr = read(fid["Components"][comp_keys[i]]["Geometry"]["rr"])
    refTag = read(fid["Components"][comp_keys[i]]["RefTag"])
    elType = read(fid["Components"][comp_keys[i]]["ElType"])
    ee_te = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["e_te"])
    rr_te = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["rr_te"])
    neigh_te = read(fid["Components"][comp_keys[i]]["Trailing_Edge"]["neigh_te"])
    n_pan = size(ee,2)
    n_vert = size(rr, 2)

    newComp = component(name, elType, refTag, ee, rr, neigh, ee_te, rr_te, neigh_te, n_pan, n_vert)
    global comps = append!(comps, [newComp])
    global nPan = nPan + n_pan
    global nVert = nVert + n_vert
end

# intitialize linear system
uinf = [30; 0; 0]
Pinf = 0
rhoinf = 1.225

# get pts in global reference frame

ee_all = zeros(4, nPan)
rr_all = zeros(3, nVert)


last_pan = 0
last_vert = 0
for i = 1:numComp
    res = findall(x -> x.==comps[i].refName, ref_keys)[1]
    orig = refs[res].origin
    orient = refs[res].orient

    # append panels
    ee_all[ :, last_pan+1:last_pan+comps[i].n_pan] = comps[i].ee .+ last_vert
    
    # rotate vertices accordingly
    rr_temp = comps[i].rr
    for j = 1:comps[i].n_vert
        rr_temp[:,j] = orient*rr_temp[:,j] + orig
    end

    rr_all[:, last_vert+1:last_vert+comps[i].n_vert] = rr_temp
    global last_pan = last_pan + comps[i].n_pan
    global last_vert = comps[i].n_vert

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


function dub(pts, loc, normal)
    cpt = mean(pts, dims=2)
    radius = norm(loc.-cpt)
    e3 = normal

    # get edge vectors
    edge_vec = zeros(3,4)
    edge_vec[:,1] = pts[:,2] .- pts[:,1]
    edge_vec[:,2] = pts[:,3] .- pts[:,2]
    edge_vec[:,3] = pts[:,4] .- pts[:,3]
    edge_vec[:,4] = pts[:,1] .- pts[:,4]

    zQ = sum((loc.-cpt).*e3)
    #display(zQ)
    dou = 0.0
    
    n_ver = 4
    for i1 = 1:n_ver
        indp = 1+i1 % n_ver
        indm1 = n_ver - (n_ver-i1+1) % n_ver
        ei = - (loc - pts[:,i1])
        ei = ei / norm(ei)
        ap1 = edge_vec[:,i1] - ei .* sum(ei .* edge_vec[:,i1])
        am1 = -edge_vec[:,indm1] .+ ei .* sum(ei .* edge_vec[:,indm1])
        ap1n = norm(ap1)
        am1n = norm(am1)
        sinB = sum(ei.*cross(am1, ap1)) / (ap1n.*am1n)
        cosB = sum(am1.*ap1 / (ap1n.*am1n))
        beta = atan(sinB, cosB)
        dou = dou + beta
    end
    #Correct the result to obtain the solid angle (from Gauss-Bonnet theorem)
    if dou < (n_ver-2) * pi + 1e-5
        dou = dou + (n_ver-2)*pi
    elseif dou > (n_ver-2) * pi - 1e-5
        dou = dou - (n_ver-2)*pi
    end
        
    return dou
end

A = zeros(nPan, nPan)
for i = 1:nPan
    for j = 1:nPan
        if i == j
            dou = -2*pi
        else
            dou = dub(rr_all[:,ee_all[:,j]], cent[:,i], normals[:,i])
        end
        A[i,j] = -dou



    end
end