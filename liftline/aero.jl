function vel_line(p1, p2, loc, gam)
    pts = [p1 p2]
    edge_vec = pts[:,2] .- pts[:,1]
    edge_len = norm(edge_vec)
    edge_uni = edge_vec ./ edge_len

    av = loc .- pts[:,1]
    ai = dot(av,edge_uni)

    R1 = norm(av)
    R2 = norm(loc.-pts[:,2])
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
    return vdou.*gam/4/pi
end

function vel_hshoe!(vel, dwash, pts, pt_idx, loc, gam)
    for i = 1:5
        v = vel_line(pts[:,pt_idx[i]], pts[:,pt_idx[i+1]], loc, gam)
        vel .= vel .+ v
        if i !=3 
            dwash .= dwash .+ v
        end
    end
end

function solve_liftline(pansys, vel_inf, rho)
    # Source: Katz, Low-Speed Aerodynamics
    rho = 1.225
    v_mag = norm(vel_inf)
    npan = size(pansys.pan_con, 2)
    pan_cpt, pan_norm, pan_area, pan_edgeuni, pan_edgelen, pan_tan = calc_liftline_props(pansys.pan_vert, pansys.pan_con)
    A = zeros(npan, npan) # influence coefficient, perpendicular to collocation point
    B = zeros(npan, npan) # influence coefficient, wake perpendicular to collocation point
    wake_inf_tan = zeros(npan, npan) # influence coefficient, wake tangent to collocation point
    RHS = zeros(npan)
    for i = 1:npan
        for j = 1:npan
            vel = zeros(3)
            dwash = zeros(3)
            vel_hshoe!(vel, dwash, pansys.pan_vert, pansys.pan_con[:,j], pan_cpt[:,i], 1)
            A[i,j] = dot(vel, pan_norm[:,i])
            B[i,j] = dot(dwash, pan_norm[:,i])
            wake_inf_tan[i,j] = dot(dwash, pan_tan[:,i])
        end
    end
    for i = 1:npan
        RHS[i] = -dot(vel_inf, pan_norm[:,i])
    end

    gam_sol = A\RHS
    gam_sol .= -gam_sol
    
    # downwash velocity from wake
    w = B*gam_sol

    # lift and drag using kutta joukouski on global flow
    dL = zeros(npan)
    dD = zeros(npan)
    for i = 1:npan
        dL[i] = rho*v_mag*gam_sol[i]*pan_edgelen[3,i]
        dD[i] = rho*w[i]*gam_sol[i]*pan_edgelen[3,i] # small angle assumption?
    end
    L = sum(dL)
    D = sum(dD)

    return L, D
end

function solve_liftline2(pansys, vel_inf, rho)
    # Modified from Katz, Low-Speed Aerodynamics
    v_mag = norm(vel_inf)
    npan = size(pansys.pan_con, 2)
    pan_cpt, pan_norm, pan_area, pan_edgeuni, pan_edgelen, pan_tan = calc_liftline_props(pansys.pan_vert, pansys.pan_con)
    A = zeros(npan, npan) # influence coefficient, perpendicular to collocation point
    wake_inf_perp = zeros(npan, npan) # influence coefficient, wake perpendicular to collocation point
    wake_inf_tan = zeros(npan, npan) # influence coefficient, wake tangent to collocation point
    RHS = zeros(npan)
    for i = 1:npan
        for j = 1:npan
            vel = zeros(3)
            wakevel = zeros(3)
            vel_hshoe!(vel, wakevel, pansys.pan_vert, pansys.pan_con[:,j], pan_cpt[:,i], 1)
            A[i,j] = dot(vel, pan_norm[:,i])
            wake_inf_perp[i,j] = dot(wakevel, pan_norm[:,i])
            wake_inf_tan[i,j] = dot(wakevel, pan_tan[:,i])
        end
    end
    for i = 1:npan
        RHS[i] = -dot(vel_inf, pan_norm[:,i])
    end
    gam_sol = A\RHS
    gam_sol .= -gam_sol

    # local tangent and perpendicular velocities
    v_perp = zeros(npan)
    v_tan = zeros(npan)
    alpha_eff = zeros(npan)
    v_total = zeros(npan)
    
    # downwash velocity from wake
    wakevel_perp = wake_inf_perp*gam_sol # perpendicular
    wakevel_tan = wake_inf_tan*gam_sol # tangent

    for i = 1:npan
        v_perp[i] = -wakevel_perp[i] + dot(vel_inf, pan_norm[:,i])
        v_tan[i] = wakevel_tan[i] + dot(vel_inf, pan_tan[:,i])
        alpha_eff[i] = atan(v_perp[i], v_tan[i])
        v_total[i] = sqrt(v_perp[i]^2 + v_tan[i]^2)
    end

    # lift and drag using kutta joukouski on local flow
    dN = zeros(npan)
    dT = zeros(npan)
    dF = zeros(3,npan)
    for i = 1:npan
        dL_local = rho*v_total[i]*gam_sol[i]*pan_edgelen[3,i]
        dN[i] = cos(alpha_eff[i])*dL_local
        dT[i] = -sin(alpha_eff[i])*dL_local
        dF[:,i] .= dN[i].*pan_norm[:,i] .+ dT[i].*pan_tan[:,i]
    end

    F = sum(dF, dims=2)

    return F
end