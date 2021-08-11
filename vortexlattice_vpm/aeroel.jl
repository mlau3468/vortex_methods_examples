using Statistics

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
    ptvel :: Array{Float64,2}
end

struct wakeRing
    gam :: Array{Float64,1}
    pts :: Array{Float64,2}
    ptvel :: Array{Float64,2}
end

struct vortRing
    gam :: Array{Float64,1}
    last_gam :: Array{Float64,1}
    dgdt :: Array{Float64,1}
    pts :: Array{Float64,2}
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
    vcpt :: Array{Float64,1} # collocatiopn point velocity
    vpts :: Array{Float64,2} # velocity of vertices
    compidx :: Array{Int64,1} # index of component panel is associated with
end

function initWakeLine(pts)
    cpt =  (pts[:,1] .+ pts[:,2])/2
    ptvel = zeros(3,2)
    gam = [0]
    return wakeLine(pts, cpt, gam, ptvel)
end

function initWakeRing(pts)
    gam = [0.0]
    ptvel = zeros(3,4)
    return wakeRing(gam, pts, ptvel)
end

function initVortRing(pts)
    gam = [0.0]
    last_gam = [0.0]
    dgdt = [0.0]
    area = [norm(cross(pts[:,2].-pts[:,1],pts[:,4].-pts[:,1]))]
    tanj_vec = pts[:,2]-pts[:,1]
    tanj_len = [norm(tanj_vec)]
    tanj_uvec = tanj_vec./tanj_len
    tani_vec = pts[:,4]-pts[:,1]
    tani_len = [norm(tani_vec)]
    tani_uvec = tani_vec./tani_len
    normal = quadNorm(pts)
    wake_vel = zeros(3) # wake induced velocity
    dp = [0.0] # pressure differential
    df = zeros(3) # force
    vcpt = zeros(3)
    vpts = zeros(3,4)
    cpt = (pts[:,1] .+ pts[:,2] .+ pts[:,3] .+ pts[:,4])./4

    compidx = [0]

    return vortRing(gam, last_gam, dgdt, pts, area, tanj_vec, tanj_len, tanj_uvec, tani_vec, tani_len, tani_uvec, normal, wake_vel, dp, df, cpt, vcpt, vpts, compidx)
end

function elemVel(panels, particles, wakelines, wakerings, loc)
    vel = [0.0;0.0;0.0]
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

function wakeElemVel(particles, wakelines, wakerings, loc)
    vel = [0.0;0.0;0.0]
    for j = 1:length(wakelines)
        vel = vel .+ velVortLine(wakelines[j], loc)
    end
    for j = 1:length(particles)
        vel = vel .+ velVortPart(particles[j], loc)
    end
    for j = 1:length(wakerings)
        vel = vel .+ velVortRing(wakerings[j], loc)
    end
    return vel

end

function vrtxline(p1, p2, p, gam)
    pts = [p1 p2]
    edge_vec = pts[:,2] .- pts[:,1]
    edge_len = norm(edge_vec)
    edge_uni = edge_vec ./ edge_len

    av = p .- pts[:,1]
    ai = dot(av,edge_uni)

    R1 = norm(av)
    R2 = norm(p.-pts[:,2])
    hv = av .- ai.*edge_uni
    hi = norm(hv)
    r_rankine = 0.001
    r_cutoff = 0.0001
    if hi > edge_len.*r_rankine
        vdou = ((edge_len.-ai)./R2 .+ai./R1) ./ (hi.^2) .* cross(edge_uni, hv)
    else
        if R1 > edge_len.*r_cutoff && R2 > edge_len.*r_cutoff
            r_ran = r_rankine .* edge_len
            vdou = ((edge_len.-ai)./R2 .+ai./R1) ./ (r_ran.^2) .* cross(edge_uni, hv)
        else
            vdou = [0.0;0.0;0.0]
        end
    end
    vel = vdou .* gam[1]
    vel = vel /4/pi
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

function stepWake!(panels, particles, wakelines, wakerings, uinf, dt)
    # calculate induced velocities at existing wake points
    for i = 1:length(wakerings)
        for j = 1:4
            vel = elemVel(panels, particles, wakelines, wakerings,wakerings[i].pts[:,j]) .+ uinf
            wakerings[i].ptvel[:,j] = vel
        end
    end

    for i = 1:length(particles)
        vel = elemVel(panels, particles, wakelines, wakerings, particles[i].cpt) .+ uinf
        particles[i].vel[:] = vel
    end

    for i = 1:length(wakelines)
        for j = 1:2
            vel = elemVel(panels, particles, wakelines, wakerings,wakelines[i].pts[:,j]) .+ uinf
            wakelines[i].ptvel[:,j] = vel
        end
    end

    # move existing wake
    for i = 1:length(wakerings)
        for j = 1:4
            wakerings[i].pts[:,j] = wakerings[i].pts[:,j] .+ wakerings[i].ptvel[:,j].*dt
        end
    end

    for i = 1:length(particles)
        particles[i].cpt[:] = particles[i].cpt .+ particles[i].vel .*dt
    end

    for i = 1:length(wakelines)
        for j = 1:2
            wakelines[i].pts[:,j] = wakelines[i].pts[:,j] .+ wakelines[i].ptvel[:,j].*dt
        end
    end
end

function shedParticles(wakeend, te_neigh, te_neighdir, wakelines)
    new_particles = []
    new_wakelines = []
    # convert wakerings into particles
    for j = 1:length(wakeend)
        
        pt1 = wakeend[j].pts[:,1]
        pt2 = wakeend[j].pts[:,2]
        pt3 = wakeend[j].pts[:,3]
        pt4 = wakeend[j].pts[:,4]
        partvec = zeros(3)

        # left side
        dir = -pt1.+pt4
        n = te_neigh[4,j]
        nd = te_neighdir[4,j]
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
        n = te_neigh[2,j]
        nd = te_neighdir[2,j]
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
        new_wakeline = initWakeLine([pt1 pt2])
        new_wakeline.gam[1] = wakeend[j].gam[1]

        push!(new_particles, new_part)
        push!(new_wakelines, new_wakeline)
    end
    return new_particles, new_wakelines
end