using LinearAlgebra
using DelimitedFiles
include("aeroel.jl")
include("geotools.jl")
include("vis.jl")

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
        panels[i].wake_vel[:] = wakeElemVel(particles, wakelines, wakerings, panels[i].cpt)
    end

    # pressure calculation
    for i = 1:length(panels)
        val = 0.0
        # i side, towards te
        if panels[i].neigh[1] > 0
            gam2 = panels[panels[i].neigh[1]].gam[1]
        else
            gam2 = 0.0
        end
        val = val .+ dot(uinf+panels[i].wake_vel, panels[i].tani_uvec).* (panels[i].gam[1]-gam2)./panels[i].tani_len

        # j side, perpendicular to te direction
        if panels[i].neigh[4] > 0
            gam2 = panels[panels[i].neigh[4]].gam[1]
        else
            gam2 = 0.0
        end
        val = val .+ dot(uinf+panels[i].wake_vel, panels[i].tanj_uvec).* (panels[i].gam[1]-gam2)./panels[i].tanj_len

        val = val .+ panels[i].dgdt[1]
        panels[i].dp[1] = -rho*val[1]
        panels[i].df[:] = -panels[i].dp*panels[i].area[1].*panels[i].normal
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

