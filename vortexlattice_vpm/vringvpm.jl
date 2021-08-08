using LinearAlgebra
using DelimitedFiles
include("geotools.jl")
include("vis.jl")


function elemVel(panels, particles, wakelines, wakerings, loc)
    vel = [0;0;0]
    for i = 1:length(panels)
        vel = vel .+ velVortRing(panels[i], loc)
    end
    for i = 1:length(particles)
        vel = vel .+ velVortPart(particles[i], loc)
    end
    for i = 1:length(wakelines)
        vel = vel .+ velVortLine(wakelines[i], loc)
    end
    for i =1:length(wakerings)
        vel = vel .+ velVortRing(wakerings[i], loc)
    end
    return vel
end

function vrtxline(p1, p2, p, gam)
    r1 = p.-p1
    r2 = p.-p2
    r0 = p2.-p1
    # check for singular conditions
    e = 1e-8
    if norm(r1) < e || norm(r2) < e || norm(cross(r1,r2))^2 < e
        vel = [0;0;0]
    else
        K = gam/(4*pi*norm(cross(r1,r2)).^2).*(dot(r0,r1)./norm(r1).-dot(r0,r2)./norm(r2))
        vel = K*cross(r1,r2)
    end
    return vel
end

function vrtxring(pts, p, gam)
    p1 = pts[:,1]
    p2 = pts[:,2]
    p3 = pts[:,3]
    p4 = pts[:,4]
    # clockwise panels
    vel1 = vrtxline(p1, p2, p, gam)
    vel2 = vrtxline(p2, p3, p, gam)
    vel3 = vrtxline(p3, p4, p, gam)
    vel4 = vrtxline(p4, p1, p, gam)
    return vel1 + vel2 + vel3 + vel4
end

function velVortRing(vring, loc)
    return vrtxring(vring.pts, loc, vring.gam[1])
end

function velVortPart(particle, loc)
    vortex_rad = 0.1
    dist = loc .- particle.cpt
    #Rosenhead kernel regularized velocity
    vvort = cross(particle.dir, dist) ./ (sqrt(sum(dist.^2)+vortex_rad.^2)).^3
    vel = vvort .* particle.gam[1]
    vel = vel /4/pi
    return vel
end

function velVortLine(vline, loc)
    vel = vrtxline(vline.pts[:,1], vline.pts[:,2], loc, vline.gam[1])
    return vel
end

function newPanGam(panel, gam, dt)
    panel.dgdt[1] = (gam - panel.gam[1])/dt
    panel.last_gam[1] = panel.gam[1]
    panel.gam[1] = gam
end

struct wakePart
    dir :: Array{Float64, 1}
    gam :: Array{Float64,1} # magnitude
    cpt :: Array{Float64,1} # position in global frame
    vel :: Array{Float64,1} # velocity in global frame
end

struct wakeLine
    pts :: Array{Float64,2}
    cpt ::Array{Float64,1}
    gam :: Array{Float64,1} # magnitude
    ptsvel :: Array{Float64,2}
end

struct wakeRing
    gam :: Array{Float64,1}
    pts :: Array{Float64,2}
    ptvel :: Array{Float64,2}
    atTE :: Array{Int64,1}
    neigh :: Array{Int64,1} # neighboring panel index
    neigh_side :: Array{Int64,1} # neighboring panel side index
    neigh_dir :: Array{Int64,1} # relative direction of neighboring segment. 1 or -1
end

struct vortRing
    gam :: Array{Float64,1}
    last_gam :: Array{Float64,1}
    dgdt :: Array{Float64,1}
    pts :: Array{Float64,2}
    vel :: Array{Float64,1}
    area :: Array{Float64,1}
    tanj_vec :: Array{Float64,1}
    tanj_len :: Array{Float64,1}
    tanj_uvec :: Array{Float64,1}
    tani_vec :: Array{Float64,1}
    tani_len :: Array{Float64,1}
    tani_uvec :: Array{Float64,1}
    normal :: Array{Float64,1}
    wake_vel :: Array{Float64,1} # wake induced velocity
    dp :: Array{Float64,1} # pressure differential
    df :: Array{Float64,1} # force
    cpt :: Array{Float64,1} # collocation point
    neigh :: Array{Int64,1} # neighboring panel index
    neigh_side :: Array{Int64,1} # neighboring panel side index
    neigh_dir :: Array{Int64,1} # relative direction of neighboring segment. 1 or -1
end

function createWakeLine(pts)
    cpt =  (pts[:,1] .+ pts[:,2])/2
    ptsvel = zeros(3,2)
    gam = [0]
    return wakeLine(pts, cpt, gam, ptsvel)
end

function createWakeRing(pts)
    gam = [0.0]
    ptvel = zeros(3,4)
    atTE = [1]
    neigh = [0,0,0,0]
    neigh_side = [0,0,0,0]
    neigh_dir = [0,0,0,0]
    return wakeRing(gam, pts, ptvel, atTE, neigh, neigh_side,neigh_dir)
end

function createVortRing(pts, vel)
    gam = [0.0]
    last_gam = [0.0]
    dgdt = [0.0]
    A = pts[:,3] .- pts[:,1]
    B = pts[:,2] .- pts[:,4]
    area = [norm(cross(pts[:,2].-pts[:,1],pts[:,4].-pts[:,1]))]
    tanj_vec = pts[:,2]-pts[:,1]
    tanj_len = [norm(tanj_vec)]
    tanj_uvec = tanj_vec./tanj_len
    tani_vec = pts[:,4]-pts[:,1]
    tani_len = [norm(tani_vec)]
    tani_uvec = tani_vec./tani_len
    normal = cross(A,B)/norm(cross(A,B))
    wake_vel = zeros(3) # wake induced velocity
    dp = [0.0] # pressure differential
    df = zeros(3) # force
    # place collocation point at 3 quarter chord
    cpt = 0.75.*(pts[:,3]+pts[:,4])/2 .+ 0.25 .* (pts[:,1]+pts[:,2])/2 # collocation point
    # place vortex ring at quarter chord, 2D kutta condition satisfid along the chord
    mcv = ((pts[:,4]-pts[:,1]) .+ (pts[:,3]-pts[:,2]))/2 # mean chord vector
    vrpts = zeros(3,4) # vortex ring points
    neigh = [0;0;0;0]
    neigh_side = [0;0;0;0]
    neigh_dir = [0;0;0;0]
    for i = 1:size(pts,2)
        vrpts[:,i] = pts[:,i] + 0.25*mcv
    end

    return vortRing(gam, last_gam, dgdt, vrpts, vel, area, tanj_vec, tanj_len, tanj_uvec, tani_vec, tani_len, tani_uvec, normal, wake_vel, dp, df, cpt, neigh, neigh_side,neigh_dir)
end

nspan = 13
nchord = 4

chord = 1
span = 8

S = span*chord
U = 50
alpha = 5

rho = 1.225

#dt = 0.1
dt = chord/U/0.1

tewidth = nspan
tsteps = 100
te_scale = 0.3
maxwakelen = 2
prefix = "test/_wing"

uinf = [U*cos(deg2rad(alpha)); 0; U*sin(deg2rad(alpha))]

# create geometry
panels = []
particles = []
wakelines = []
wakerings = []
te_idx = []
wakelen = 0
for i = 0:nchord-1
    for j = 0:nspan-1
        p1 = [i*chord/nchord; j*span/nspan; 0]
        p2 = [i*chord/nchord; (j+1)*span/nspan; 0]
        p3 = [(i+1)*chord/nchord; (j+1)*span/nspan; 0]
        p4 = [(i+1)*chord/nchord; j*span/nspan; 0]
        pts = [p1 p2 p3 p4]
        vel = [0, 0, 0]
        new_pan = createVortRing(pts, vel)
        push!(panels, new_pan)
        if i+1==nchord
            push!(te_idx, i*nspan + j + 1)
        end
    end
   
end

getNeighbors!(panels)

A = zeros(length(panels), length(panels))
RHS = zeros(length(panels))

for i = 1:length(panels)
    for j = 1:length(panels)
        # influence of jth panel on ith collocation point
        vel = vrtxring(panels[j].pts, panels[i].cpt, 1)
        A[i,j] = dot(vel, panels[i].normal)
    end
end

for t = 1:tsteps

    # build rhs vector
    RHS[:] .= 0.0
    for i = 1:length(panels)
        RHS[i] = -dot(uinf, panels[i].normal)
        RHS[i] = RHS[i] - dot(panels[i].wake_vel, panels[i].normal)
    end

    # solve matrix for panel gamma
    sol = A\RHS
    for i = 1:length(panels)
        newPanGam(panels[i], sol[i], dt)
    end

    # calculate new wake elements
    new_particles = []
    new_wakelines = []
    new_wakerings = []
    for i = 1:length(te_idx)
        idx = te_idx[i]
        p1 = panels[idx].pts[:,4]
        p2 = panels[idx].pts[:,3]

        vel1 = elemVel(panels, particles, wakelines, wakerings, p1) .+ uinf
        vel2 = elemVel(panels, particles, wakelines, wakerings, p2) .+ uinf
        gam = panels[idx].last_gam[1]
        p3 = p2 .+ te_scale.*vel2.*dt
        p4 = p1 .+ te_scale.*vel1.*dt
        new_wakering = createWakeRing([p1 p2 p3 p4])
        new_wakering.gam[1] = gam
        push!(new_wakerings, new_wakering)
    end

    # calculate induced velocities at existing wake points
    for i = 1:length(wakerings)
        for j = 1:4
            vel = elemVel(panels, particles, wakelines, wakerings, wakerings[i].pts[:,j]) .+ uinf
            wakerings[i].ptvel[:,j] = vel
        end
    end

    for i = 1:length(particles)
        vel = elemVel(panels, particles, wakelines, wakerings, particles[i].cpt) .+ uinf
        particles[i].vel[:] = vel
    end

    # move existing wake
    for i = 1:length(wakerings)
        if wakerings[i].atTE[1] == 1
            for j = 1:2
            wakerings[i].pts[:,j] = wakerings[i].pts[:,j] .+ wakerings[i].ptvel[:,j].*dt.*te_scale
            end
            for j = 3:4
                wakerings[i].pts[:,j] = wakerings[i].pts[:,j] .+ wakerings[i].ptvel[:,j].*dt
                end
            wakerings[i].atTE[1] = 0
        else
            for j = 1:4
                wakerings[i].pts[:,j] = wakerings[i].pts[:,j] .+ wakerings[i].ptvel[:,j].*dt
            end
        end
    end

    for i = 1:length(particles)
        particles[i].cpt[:] = particles[i].cpt .+ particles[i].vel .*dt
    end

    # add new wakerings
    global wakerings = cat(new_wakerings, wakerings, dims=1)

    global wakelen = wakelen + 1
     # convert end wake rings into particles
    if wakelen > maxwakelen
        wakeend = wakerings[tewidth*(wakelen-1)+1:end]
        for j = 1:length(wakeend)
            getNeighbors!(wakeend)

            pt1 = wakeend[j].pts[:,1]
            pt2 = wakeend[j].pts[:,2]
            pt3 = wakeend[j].pts[:,3]
            pt4 = wakeend[j].pts[:,4]
            partvec = zeros(3)

            # left side
            dir = -pt1.+pt4
            n = wakeend[j].neigh[4]
            nd = wakeend[j].neigh_dir[4]
            if n > 0
                ave = wakeend[j].gam[1] -  nd*wakeend[n].gam[1]
                ave = ave/2
            else
                ave = wakeend[j].gam[1]
            end
            partvec = partvec .+ dir.*ave

            # right side
            dir = -pt3 .+ pt2
            # if has neighboring panel
            n = wakeend[j].neigh[2]
            nd = wakeend[j].neigh_dir[2]
            if n > 0
                ave = wakeend[j].gam[1] -  nd*wakeend[n].gam[1]
                ave = ave/2
            else
                ave = wakeend[j].gam[1]
            end
            partvec = partvec .+ dir.*ave
    
            # end side
            dir = -pt4 .+ pt3
            if length(wakelines) > 0
                ave = wakeend[j].gam[1] - wakelines[j].gam[1]
            else
                ave = wakeend[j].gam[1] - 0
            end
            partvec = partvec .+ dir.*ave
    
            # create particle
            posp = (pt1 .+ pt2 .+ pt3 .+ pt4)./4
            mag_partvec = norm(partvec)

            if wakeend[j].gam[1] < 1e-13
                dir = partvec
            else
                dir = partvec./mag_partvec
            end
            vel = [0;0;0] # particle released, floating.
            new_part = wakePart(dir, [mag_partvec], posp, vel)
    
            # new wake line left over from wake panel to particle conversion.
            # inject downstream for the next iteration
            new_wakeline = createWakeLine([pt4 pt3])
            new_wakeline.gam[1] = wakeend[j].gam[1]
    
            push!(new_particles, new_part)
            push!(new_wakelines, new_wakeline)
        end
        global wakelen = wakelen - 1
        global wakerings = wakerings[1:tewidth*wakelen]

        # append wakelines to system 
        if length(wakelines) == tewidth
            for k = 1:length(wakelines)
                wakelines[k] = new_wakelines[k]
            end
        else
            global wakelines = new_wakelines
        end

        # append particles to system
        global particles = cat(particles, new_particles, dims=1)
        
    end

    panels2vtk(panels, prefix * "_panels_$t.vtu")
    particles2vtk(particles, prefix * "_particles_$t.vtu")
    wakepanels2vtk(wakerings, prefix * "_wakerings_$t.vtu")

    # update panel wake_vel
    for i =1:length(panels)
        vel = [0;0;0]
        for j = 1:length(wakelines)
            vel = vel .+ velVortLine(wakelines[j], panels[i].cpt)
        end
        for j = 1:length(particles)
            vel = vel .+ velVortPart(particles[j], panels[i].cpt)
        end
        for j = 1:length(wakerings)
            vel = vel .+ velVortRing(wakerings[j], panels[i].cpt)
        end
        panels[i].wake_vel[:] = vel
    end

    # pressure calculation
    for i = 0:nchord-1
        for j = 0:nspan-1
            idx = i*nspan + j + 1
            val = 0.0
            if i > 0
                #idx2 = idx - nspan
                idx2 = (i-1)*nspan+j+1
                gam2 = panels[idx2].gam[1]
            else
                gam2 = 0.0
            end
            val = val .+ dot(uinf+panels[idx].wake_vel, panels[idx].tani_uvec).* (panels[idx].gam[1]-gam2)./panels[idx].tani_len

            if j > 0
                #idx2 = idx - nchord
                idx2 = i*nspan+(j-1)+1
                gam2 = panels[idx2].gam[1]
            else
                gam2 = 0.0
            end
            val = val .+ dot(uinf+panels[idx].wake_vel, panels[idx].tanj_uvec).* (panels[idx].gam[1]-gam2)./panels[idx].tanj_len

            val = val .+ panels[idx].dgdt[1]
            panels[idx].dp[1] = -rho*val[1]
            panels[idx].df[:] = -panels[idx].dp*panels[idx].area[1].*panels[idx].normal            
        end
    end

    # total forces
    total_force = [0;0;0]
    for i=1:length(panels)
        total_force = total_force .+ panels[i].df
    end
    # lift coefficient
    cl = cos(deg2rad(alpha))*total_force[3] - sin(deg2rad(alpha)) * total_force[1]
    cl = cl/(1/2*rho*U^2*S)
    println("Step: $t, CL=$cl")


end

