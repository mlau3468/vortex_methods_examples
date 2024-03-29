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
    r_rankine = 0.0001
    r_cutoff = 0.00001
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

function solve_liftline_katz(pansys, vel_inf, rho)
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
    
    # downwash velocity from wake
    w = .-B*gam_sol

    # lift and drag using kutta joukouski on global flow
    dL = zeros(npan)
    dD = zeros(npan)
    for i = 1:npan
        dL[i] = -rho*v_mag*gam_sol[i]*pan_edgelen[3,i]
        dD[i] = rho*w[i]*gam_sol[i]*pan_edgelen[3,i] # small angle assumption?
    end
    L = sum(dL)
    D = sum(dD)

    return L, D
end

function solve_liftline_katz_mod(pansys, vel_inf, rho)
    # Modified from Katz, Low-Speed Aerodynamics
    npan = size(pansys.pan_con, 2)
    pan_cpt, pan_norm, pan_area, pan_edgeuni, pan_edgelen, pan_tan = calc_liftline_props(pansys.pan_vert, pansys.pan_con)
    A = zeros(npan, npan) # influence coefficient, perpendicular to collocation point
    B = zeros(npan, npan) # influence coefficient, tangent to collocation point
    wake_inf_perp = zeros(npan, npan) # influence coefficient, wake perpendicular to collocation point
    wake_inf_tan = zeros(npan, npan) # influence coefficient, wake tangent to collocation point
    RHS = zeros(npan)

    # local tangent and perpendicular velocities
    v_perp = zeros(npan)
    v_tan = zeros(npan)
    alpha_eff = zeros(npan)
    alpha_ind = zeros(npan)
    v_total = zeros(npan)
    dwash = zeros(npan)

    for i = 1:npan
        for j = 1:npan
            vel = zeros(3)
            wakevel = zeros(3)
            vel_hshoe!(vel, wakevel, pansys.pan_vert, pansys.pan_con[:,j], pan_cpt[:,i], 1)
            A[i,j] = dot(vel, pan_norm[:,i])
            B[i,j] = dot(vel, pan_tan[:,i])
            wake_inf_perp[i,j] = dot(wakevel, pan_norm[:,i])
            wake_inf_tan[i,j] = dot(wakevel, pan_tan[:,i])
        end
    end
    for i = 1:npan
        RHS[i] = -dot(vel_inf, pan_norm[:,i])
    end
    gam_sol = A\RHS
    
    # calculate effective angle of attack
    for i = 1:npan
        v_perp[i] = sum(wake_inf_perp[i,:].*gam_sol) + dot(vel_inf, pan_norm[:,i]) # wake induced + freestream
        v_tan[i] = sum(wake_inf_tan[i,:].*gam_sol) + dot(vel_inf, pan_tan[:,i]) # wake induced + freestream
        v_total[i] = norm(vel_inf)
        dwash[i] = sum(wake_inf_perp[i,:].*gam_sol)
        alpha_eff[i] = atan(v_perp[i], v_tan[i])
        alpha_ind[i] = atan(-dwash[i], v_total[i])
    end

    # lift and drag using kutta joukouski on local flow
    dN = zeros(npan)
    dT = zeros(npan)
    dF = zeros(3,npan)
    cls = zeros(npan)
    for i = 1:npan
        dL_local = -rho*v_total[i]*gam_sol[i]*pan_edgelen[3,i]
        #cls[i] = -2*gam_sol[i]/
        #cl_local = 2*pi*alpha_eff[i]
        #dL_local = 1/2*rho*v_total[i]^2*pan_area[i]*cl_local
        dN[i] = cos(alpha_eff[i])*dL_local
        dT[i] = -sin(alpha_eff[i])*dL_local
        dF[:,i] .= dN[i].*pan_norm[:,i] .+ dT[i].*pan_tan[:,i]
    end

    F = sum(dF, dims=2)

    return F
end

function solve_liftline_weissinger(pansys, vel_inf, rho, affile)
    # Source: Owens, Weissinger's Modle of the Nonlinear Lifting-Line MEthod for Aircraft Design
    npan = size(pansys.pan_con, 2)
    pan_cpt, pan_norm, pan_area, pan_edgeuni, pan_edgelen, pan_tan = calc_liftline_props(pansys.pan_vert, pansys.pan_con)
    A = zeros(npan, npan) # influence coefficient, perpendicular to collocation point
    wake_inf_perp = zeros(npan, npan) # influence coefficient, wake perpendicular to collocation point
    wake_inf_tan = zeros(npan, npan) # influence coefficient, wake tangent to collocation point
    RHS = zeros(npan)

    # local tangent and perpendicular velocities
    v_perp = zeros(npan)
    v_tan = zeros(npan)
    alpha_eff = zeros(npan)
    v_total = zeros(npan)
    cls = zeros(npan)
    cds = zeros(npan)
    residuals = zeros(npan)
    tol = 1e-7
    max_iter = 2000
    iter = 1
    rlx = 0.01

    # airfoil polar
    af_data = readdlm(affile, ',')
    cl_interp = LinearInterpolation(af_data[:,1], af_data[:,2])
    cd_interp = LinearInterpolation(af_data[:,1], af_data[:,3])

    # Inflence matrices
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

    # Initial solution for circulation, uncorrected
    gam_sol = A\RHS
    # Initial solution for ciculation at 0
    #gam_sol = zeros(npan)
    gam_temp = zeros(npan)

    done = false
    while !done
        # calculate effective angle of attack, cl, cd
        for i = 1:npan
            v_perp[i] = sum(wake_inf_perp[i,:].*gam_sol) + dot(vel_inf, pan_norm[:,i]) # wake induced + freestream
            v_tan[i] = sum(wake_inf_tan[i,:].*gam_sol) + dot(vel_inf, pan_tan[:,i]) # wake induced + freestream
            v_total[i] = sqrt(v_perp[i]^2 + v_tan[i]^2)
            v_total[i] = norm(vel_inf)

            chord = (pan_edgelen[2,i] + pan_edgelen[4,i])/2
            
            alpha_eff[i] = atan(v_perp[i], v_tan[i])
            cls[i] = cl_interp(rad2deg(alpha_eff[i]))
            cds[i] = cd_interp(rad2deg(alpha_eff[i]))
            gamnew = -0.5*v_total[i]*cls[i]*chord
            residuals[i] = abs(gamnew - gam_sol[i])
            gam_temp[i] = gam_sol[i]*(1-rlx) + rlx*gamnew
        end

        iter += 1

        if maximum(residuals) < tol
            done = true
        elseif iter == max_iter
            done = true
            println("Maximum iterations exceeded.")
        else
            gam_sol[:] .= gam_temp
        end
    end
    
    # lift and drag using kutta joukouski on local flow
    dN = zeros(npan)
    dT = zeros(npan)
    dF = zeros(3,npan)
    for i = 1:npan
        dL_local = -rho*v_total[i]*gam_sol[i]*pan_edgelen[3,i]
        dN[i] = cos(alpha_eff[i])*dL_local
        dT[i] = -sin(alpha_eff[i])*dL_local
        dF[:,i] .= dN[i].*pan_norm[:,i] .+ dT[i].*pan_tan[:,i]
    end
    F_inv = sum(dF, dims=2)

    # lift and drag using airfoil polars
    for i = 1:npan
        dL_local = 0.5*rho*v_total[i]^2*pan_area[i]*cls[i]
        dD_local = 0.5*rho*v_total[i]^2*pan_area[i]*cds[i]
        dN[i] = cos(alpha_eff[i])*dL_local + sin(alpha_eff[i])*dD_local
        dT[i] = -sin(alpha_eff[i])*dL_local + cos(alpha_eff[i])*dD_local
        dF[:,i] .= dN[i].*pan_norm[:,i] .+ dT[i].*pan_tan[:,i]
    end
    F_visc = sum(dF, dims=2)
    return F_inv, F_visc
end


function solve_liftline_vandam(pansys, vel_inf, rho, affile)
    # Source: van Dam, The aerodynamic design of multi-element high-lift systems for transport airplanes
    npan = size(pansys.pan_con, 2)
    pan_cpt, pan_norm, pan_area, pan_edgeuni, pan_edgelen, pan_tan = calc_liftline_props(pansys.pan_vert, pansys.pan_con)
    A = zeros(npan, npan) # influence coefficient, perpendicular to collocation point
    wake_inf_perp = zeros(npan, npan) # influence coefficient, wake perpendicular to collocation point
    wake_inf_tan = zeros(npan, npan) # influence coefficient, wake tangent to collocation point
    RHS = zeros(npan)

    # local tangent and perpendicular velocities
    v_perp = zeros(npan)
    v_tan = zeros(npan)
    alpha_eff = zeros(npan)
    v_total = zeros(npan)
    residuals = zeros(npan)
    tol = 1e-6
    max_iter = 500
    iter = 1
    rlx = 1

    # airfoil polar
    af_data = readdlm(affile, ',')
    cl_interp = LinearInterpolation(af_data[:,1], af_data[:,2])
    cd_interp = LinearInterpolation(af_data[:,1], af_data[:,3])

    vel_mag = norm(vel_inf)
    alphas = zeros(npan)
    dalpha = zeros(npan)
    cl_inv = zeros(npan)
    cl_visc = zeros(npan)
    cd_visc = zeros(npan)
    gam_sol = zeros(npan)

    # Local angle of attack
    for i = 1:npan
        v_p = dot(vel_inf, pan_norm[:,i]) # freestream
        v_t = dot(vel_inf, pan_tan[:,i]) # freestream
        alphas[i] = atan(v_p, v_t)
    end

    # Inflence matrices
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

    done = false
    while !done
        
        for i = 1:npan
            RHS[i] = -sin(alphas[i] - dalpha[i])*vel_mag
        end
    
        gam_sol[:] = A\RHS
    
        # calculate effective angle of attack, cl, cd
        for i = 1:npan
            v_perp[i] = sum(wake_inf_perp[i,:].*gam_sol) + dot(vel_inf, pan_norm[:,i]) # wake induced + freestream
            v_tan[i] = sum(wake_inf_tan[i,:].*gam_sol) + dot(vel_inf, pan_tan[:,i]) # wake induced + freestream
            v_total[i] = norm(vel_inf)

            chord = (pan_edgelen[2,i] + pan_edgelen[4,i])/2
            
            alpha_eff[i] = atan(v_perp[i], v_tan[i])
            cl_inv[i] = -2*gam_sol[i]/v_total[i]/chord
            cl_visc[i] = cl_interp(rad2deg(alpha_eff[i]))
            cd_visc[i] = cd_interp(rad2deg(alpha_eff[i]))
            dalpha_new = (cl_inv[i] - cl_visc[i])/(2*pi)
            dalpha[i] = dalpha[i] + rlx*dalpha_new
            residuals[i] = abs(cl_visc[i] - cl_inv[i])
        end

        iter += 1
        if maximum(residuals) < tol
            done = true
        elseif iter == max_iter
            done = true
            println("Maximum iterations exceeded.")
        end
    end

    # lift and drag using kutta joukouski on local flow
    dN = zeros(npan)
    dT = zeros(npan)
    dF = zeros(3,npan)
    for i = 1:npan
        dL_local = -rho*v_total[i]*gam_sol[i]*pan_edgelen[3,i]
        dN[i] = cos(alpha_eff[i])*dL_local
        dT[i] = -sin(alpha_eff[i])*dL_local
        dF[:,i] .= dN[i].*pan_norm[:,i] .+ dT[i].*pan_tan[:,i]
    end
    F_inv = sum(dF, dims=2)

    # lift and drag using airfoil polars
    for i = 1:npan
        dL_local = 0.5*rho*v_total[i]^2*pan_area[i]*cl_visc[i]
        dD_local = 0.5*rho*v_total[i]^2*pan_area[i]*cd_visc[i]
        dN[i] = cos(alpha_eff[i])*dL_local + sin(alpha_eff[i])*dD_local
        dT[i] = -sin(alpha_eff[i])*dL_local + cos(alpha_eff[i])*dD_local
        dF[:,i] .= dN[i].*pan_norm[:,i] .+ dT[i].*pan_tan[:,i]
    end
    F_visc = sum(dF, dims=2)
    return F_inv, F_visc
end