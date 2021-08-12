using Interpolations

function samePt(pt1, pt2)
    tol = 1e-6
    return norm(pt1.-pt2) < tol
end

function calcneighbors(panels)
    idx = [1 2 3 4 1]
    neigh_idx = zeros(Int64,4,length(panels))
    neigh_side = zeros(Int64,4,length(panels))
    neigh_dir = zeros(Int64,4,length(panels))
    for i = 1:length(panels)
        pan1 = panels[i]
        for ii=1:4
            if neigh_idx[ii,i] == 0 # if currently no neighbor, check
                p1 = pan1.pts[:,idx[ii]]
                p2 = pan1.pts[:,idx[ii+1]]
                for j = 1:length(panels)
                    if i != j
                        pan2 = panels[j]
                        for jj = 1:4
                            p11 = pan2.pts[:,idx[jj]]
                            p22 = pan2.pts[:,idx[jj+1]]
                            if (samePt(p1,p22) && samePt(p2,p11))
                                # is neighbor
                                neigh_idx[ii,i] = j
                                neigh_idx[jj,j] = i
                                neigh_side[ii,i] = jj
                                neigh_side[jj,j] = ii
                                neigh_dir[ii,i] = -1
                                neigh_dir[jj,j] = -1
                            elseif (samePt(p1,p11) && samePt(p2,p22))
                                # is neighbor
                                neigh_idx[ii,i] = j
                                neigh_idx[jj,j] = i
                                neigh_side[ii,i] = jj
                                neigh_side[jj,j] = ii
                                neigh_dir[ii,i] = 1
                                neigh_dir[jj,j] = 1
                            end
                        end
                    end
                end
            end
        end
    end
    return neigh_idx, neigh_side, neigh_dir
end

function stepRotMat(omega_body, dt)
    rotation = zeros(Float64, 3,3)
    if norm(omega_body) > 0
        #source: https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas/rotation-about-an-arbitrary-axis-in-3-dimensions
        theta = norm(omega_body) .* dt # radians
        rvec = omega_body./norm(omega_body)
        u = rvec[1]
        v = rvec[2]
        w = rvec[3]
        a11 = u^2 + (v^2+w^2)*cos(theta)
        a12 = u*v*(1-cos(theta))-w*sin(theta)
        a13 = u*w*(1-cos(theta))+v*sin(theta)
        a21 = u*v*(1-cos(theta))+w*sin(theta)
        a22 = v^2 + (u^2+w^2)*cos(theta)
        a23 = v*w*(1-cos(theta))-u*sin(theta)
        a31 = u*w*(1-cos(theta))-v*sin(theta)
        a32 = v*w*(1-cos(theta))+u*sin(theta)
        a33 = w^2 + (u^2+v^2)*cos(theta)

        rotation = [a11 a12 a13; a21 a22 a23; a31 a32 a33]
    else
        rotation = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    end
    return rotation
end

function createRect(name, span, chord, nspan, nchord, offset, twist)
    twist = deg2rad(twist)
    twist_loc = [0.25*chord; 0; 0]
    rot = [cos(twist) 0 sin(twist); 0 1 0; -sin(twist) 0 cos(twist)]
    # create geometry
    te_idx = Int64[]
    panels = []
    for i = 0:nchord-1
        for j = 0:nspan-1
            p1 = [i*chord/nchord; j*span/nspan; 0]
            p2 = [i*chord/nchord; (j+1)*span/nspan; 0]
            p3 = [(i+1)*chord/nchord; (j+1)*span/nspan; 0]
            p4 = [(i+1)*chord/nchord; j*span/nspan; 0]
            pts = [p1 p2 p3 p4]
            for k=1:4
                p = pts[:,k]
                p = p .- twist_loc
                p = rot * p
                p = p + twist_loc
                pts[:,k] = p .+ offset
            end
            new_pan = initVortRing(pts)
            new_pan.compidx[1] = 1
            push!(panels, new_pan)
            if i+1==nchord
                push!(te_idx, i*nspan + j + 1)
            end
        end
    end
    panidx = 1:length(panels)
    newcomp = component(name, panidx, te_idx, [0;0;0], [0;0;0], [0;0;0])
    return newcomp, panels
end

struct component
    name :: String
    panidx :: Array{Int64,1} # indices of panels in global list
    teidx :: Array{Int64,1} # indices of trailing edge panels in global list
    origin :: Array{Float64,1} # origin wrt intertial
    omega :: Array{Float64,1} # rotational velocity wrt inertial
    vel :: Array{Float64,1} # velocity of origin wrt inertial
end

function quadNorm(pts)
    A = pts[:,3] .- pts[:,1]
    B = pts[:,2] .- pts[:,4]
    normal = cross(A,B)/norm(cross(A,B))
    return normal
end

function setPointVels!(components, panels)
    # set panel point velocities according to compnent rotation and linear velocity definitions
    for c = 1:length(components)
        comp = components[c]
        ombody = comp.omega
        vbody = comp.vel
        for i = 1:length(comp.panidx)
            idx = comp.panidx[i]
            for j = 1:4
                panels[idx].vpts[:,j] = cross(ombody, panels[idx].pts[:,j]) .+ vbody
            end
            panels[idx].vcpt[:] = cross(ombody, panels[idx].cpt[:]) .+ vbody
        end
    end
end

function stepGeometry!(components, panels, dt)
    for c = 1:length(components)
        comp = components[c]
        ombody = comp.omega
        vbody = comp.vel
        rotm = stepRotMat(ombody, dt)
        # update component origin
        comp.origin[:] = comp.origin .+ vbody .* dt

        # update panels of this component
        for i = 1:length(comp.panidx)
            idx = comp.panidx[i]
            # update point positions
            for j =1:4
                panels[idx].pts[:,j] = rotm*(panels[idx].pts[:,j] .+ vbody.*dt.-comp.origin) .+ comp.origin
            end
            # update collocation point and other vectors
            panels[idx].cpt[:] = (panels[idx].pts[:,1] .+ panels[idx].pts[:,2] .+ panels[idx].pts[:,3] .+ panels[idx].pts[:,4])./4
            panels[idx].normal[:] = quadNorm(panels[idx].pts)

            pts = panels[idx].pts
            tanj_vec = pts[:,2]-pts[:,1]
            tanj_len = [norm(tanj_vec)]
            tanj_uvec = tanj_vec./tanj_len
            tani_vec = pts[:,4]-pts[:,1]
            tani_len = [norm(tani_vec)]
            tani_uvec = tani_vec./tani_len
    
            panels[idx].tanj_vec[:] = tanj_vec
            panels[idx].tanj_len[:] = tanj_len
            panels[idx].tanj_uvec[:] = tanj_uvec
            panels[idx].tani_vec[:] = tani_vec
            panels[idx].tani_len[:] = tani_len
            panels[idx].tani_uvec[:] = tani_uvec

             # update point velocities
            for j = 1:4
                panels[idx].vpts[:,j] = cross(ombody, panels[idx].pts[:,j].-comp.origin) .+ vbody
            end
            panels[idx].vcpt[:] = cross(ombody, panels[idx].cpt[:].-comp.origin) .+ vbody
        end
    end
end

function test_geo(componentsin, panelsin, dt, prefix)
    panels = deepcopy(panelsin)
    components = deepcopy(componentsin)

    # set initial point velocities
    setPointVels!(components, panels)

    panels2vtk(panels, prefix * "_test_panels_0.vtu")

    for t = 1:tsteps
        stepGeometry!(components, panels, dt)

        panels2vtk(panels, prefix * "_test_panels_$t.vtu")
    end
end

function createWing(name, ys, chords, sweeps, twists, nspan, nchord, xref)
    spanspace = "cosOB"
    # convert sweep to x location
    ys2 = cat([0.0;],ys, dims=1)
    dxs = diff(ys2[1:end]).*tan.(deg2rad.(sweeps))
    xs = cumsum(dxs)
    # resample
    # create linear interpolations
    chords_itp = LinearInterpolation(ys, chords)
    sweeps_itp = LinearInterpolation(ys, sweeps)
    twists_itp = LinearInterpolation(ys, twists)
    xs_int = LinearInterpolation(ys, xs)
    
    w = 0.5
    k = sin.(LinRange(0,pi/2,nspan+1))
    k2 = LinRange(0,1,nspan+1)
    b = k.*w .+ k2.*(1-w)

    ys = ys[1] .+ (ys[end]-ys[1]).*b
    chords = chords_itp(ys)
    sweeps = sweeps_itp(ys)
    twists = twists_itp(ys)
    xs = xs_int(ys)
    
    # calculate geometry
    points = zeros(3,(nchord+1)*(nspan+1))
    faces = zeros(Int64, 4,nchord*nspan)

    panels = vortRing[]
    
    b = 1
    for j = 1:nspan+1
        xc = LinRange(-xref*chords[j]+xs[j],(1-xref)*chords[j]+xs[j],nchord+1)
        a = deg2rad(twists[j])
        rotm = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)]
        for i =1:nchord+1
            p = [xc[i];ys[j];0]
            p = rotm*p
            points[:,b] = p
            b = b + 1
        end
    end

    b = 1
    for i =1:nchord
        for j = 1:nspan
            faces[:,b] = [(j-1)*(nchord+1)+i;(j)*(nchord+1)+i;(j)*(nchord+1)+i+1;(j-1)*(nchord+1)+i+1]
            b = b + 1
        end
    end
    
    b = 1
    te_idx = Int64[]
    for i =1:nchord
        for j = 1:nspan
            pts = points[:,faces[:,b]]
            new_pan = initVortRing(pts)
            new_pan.compidx[1] = 1
            push!(panels, new_pan)
            if i == nchord
                push!(te_idx, b)
            end
            b = b + 1
        end
    end
    panidx = 1:length(panels)
    newcomp = component(name, panidx, te_idx, [0;0;0], [0;0;0], [0;0;0])
    return newcomp, panels
end