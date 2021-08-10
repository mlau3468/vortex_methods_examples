
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

function createRect(span, chord, nspan, nchord, origin, twist)
    twist = deg2rad(twist)
    twist_loc = [0.25*chord; 0; 0]
    rot = [cos(twist) 0 sin(twist); 0 1 0; -sin(twist) 0 cos(twist)]
    # create geometry
    te_idx = []
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
                pts[:,k] = p .+ origin
            end
            vel = [0; 0; 0]
            new_pan = initVortRing(pts, vel)
            push!(panels, new_pan)
            if i+1==nchord
                push!(te_idx, i*nspan + j + 1)
            end
        end
    end

    return panels, te_idx
end

function quadNorm(pts)
    A = pts[:,3] .- pts[:,1]
    B = pts[:,2] .- pts[:,4]
    normal = cross(A,B)/norm(cross(A,B))
    return normal
end

function test_geo(panelsin, origin, vbody, ombody)
    rotm = stepRotMat(ombody, dt)
    panels = deepcopy(panelsin)

    # calculate intial geometry velocities
    for i = 1:length(panels)
        for j = 1:4
            panels[i].vpts[:,j] = cross(ombody, panels[i].pts[:,j]) .+ vbody
        end
        panels[i].vcpt[:] = cross(ombody, panels[i].cpt[:]) .+ vbody
    end

    panels2vtk(panels, prefix * "_test_panels_0.vtu")

    for t = 1:tsteps
        # move geometry
        origin = origin .+ vbody.*dt
        for i = 1:length(panels)
            # update point positions
            for j =1:4
                panels[i].pts[:,j] = rotm*(panels[i].pts[:,j] .+ vbody.*dt.-origin) .+ origin
            end

            panels[i].cpt[:] = (panels[i].pts[:,1] .+ panels[i].pts[:,2] .+ panels[i].pts[:,3] .+ panels[i].pts[:,4])./4
            panels[i].normal[:] = quadNorm(panels[i].pts)

            pts = panels[i].pts

            tanj_vec = pts[:,2]-pts[:,1]
            tanj_len = [norm(tanj_vec)]
            tanj_uvec = tanj_vec./tanj_len
            tani_vec = pts[:,4]-pts[:,1]
            tani_len = [norm(tani_vec)]
            tani_uvec = tani_vec./tani_len
    
            panels[i].tanj_vec[:] = tanj_vec
            panels[i].tanj_len[:] = tanj_len
            panels[i].tanj_uvec[:] = tanj_uvec
            panels[i].tani_vec[:] = tani_vec
            panels[i].tani_len[:] = tani_len
            panels[i].tani_uvec[:] = tani_uvec

            # update point velocities
            for j = 1:4
                panels[i].vpts[:,j] = cross(panels[i].pts[:,j], ombody) .+ vbody
            end
            panels[i].vcpt[:] = cross(panels[i].cpt[:], ombody) .+ vbody
        end

        panels2vtk(panels, prefix * "_test_panels_$t.vtu")
    end
end