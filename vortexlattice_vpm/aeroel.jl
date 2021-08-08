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
    return wakeRing(gam, pts, ptvel, atTE)
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
    for i = 1:size(pts,2)
        vrpts[:,i] = pts[:,i] + 0.25*mcv
    end

    return vortRing(gam, last_gam, dgdt, vrpts, vel, area, tanj_vec, tanj_len, tanj_uvec, tani_vec, tani_len, tani_uvec, normal, wake_vel, dp, df, cpt)
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