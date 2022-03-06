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
    return vdou .* gam/4/pi
end

function vrtxring(pts, p, gam)
    p1 = pts[:,1]
    p2 = pts[:,2]
    p3 = pts[:,3]
    p4 = pts[:,4]
    # clockwise panels
    vel = zeros(3)
    vel = vel + vrtxline(p1, p2, p, gam)
    vel = vel + vrtxline(p2, p3, p, gam)
    vel = vel + vrtxline(p3, p4, p, gam)
    #vel = vel + vrtxline(p4, p1, p, gam)
    return vel
end

function hshoe(pts, p)
    # calculates influence from one horshoe vortex on point p
    # pts A,B,C,D
    # vortexlines: A-B, B-C, C-D 
    p1 = pts[:,1]
    p2 = pts[:,2]
    p3 = pts[:,3]
    p4 = pts[:,4]
    
    uvw = zeros(3)
    uvw = uvw + vrtxline(p1, p2, p, 1)
    uvw = uvw + vrtxline(p2, p3, p, 1)
    uvw = uvw + vrtxline(p3, p4, p, 1)

    dwash = zeros(3)
    dwash = dwash + vrtxline(p1, p2, p, 1)
    dwash = dwash + vrtxline(p3, p4, p, 1)
    return uvw, dwash
end