include("panel.jl")

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
newRef = refFrame("test", "0", [], false, false, [0;0;0], [0.9961947 0.0 -0.0871557; 0.0 1.0 0.0; 0.0871557 0.0 0.9961947]')
refs = append!(refs, [newRef])
ref_keys = ["test"]

# Build component lists
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

# get pts in global reference frame and get full system

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




writedlm("A_pre.csv", A)
A = factorize2(A)


writedlm("A.csv", A)
writedlm("B.csv", B)



println("finished")

#=
C = readdlm("test.txt", ',')
D = zeros(nPan, nPan)
for i=1:size(C,1)
    D[Int(C[i,1]), Int(C[i,2])] = C[i,3]
end
=#