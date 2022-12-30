function airfoil_doublet_dirichlet(pan_vert::Matrix{<:Real}, aoa::Real; compute_full_potential::Bool=false)
    # set compute_full_potential to true to compute full potential and take derivative to get velocity

    npan = size(pan_vert,2) - 1

    # Compute panel properties
    pan_cpt = calc_panel_2d_collocation_points(pan_vert) # collocation points
    pan_norm = calc_panel_2d_normals(pan_vert) # panel normals
    pan_tan = calc_panel_2d_tangents(pan_vert) # panel tangents
    pan_len = calc_panel_2d_lengths(pan_vert)

    # Freestream
    U = 1
    alf = deg2rad(aoa)
    chord = 1
    rho = 1
    pstatic = 0
    p0 = pstatic + 0.5*rho*U^2 # stagnation pressure

    u_vec = zeros(2)
    u_vec[1] = cos(alf)
    u_vec[2] = sin(alf)

    A = zeros(npan, npan) # Influence coefficent for potential just inside airfoil
    A_pot_out = zeros(npan, npan) # Influence coefficeint for potential just outside airfoil
    RHS = zeros(npan)

    # Influence coefficients of doublet
    for i = 1:npan
        for j = 1:npan
            if i == j
                A[i,j] = pot_line_doublet_2d_self() # limit approaching from inside airfoil
                A_pot_out[i,j] = pot_line_doublet_2d_self(false) # limit approaching from outside airfoil
            else
                A[i,j] = pot_line_doublet_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
                A_pot_out[i,j] = A[i,j]
            end
        end
        # trailing edge influence
        te = pot_line_doublet_2d(pan_vert[:,1], [1e3;0.0], pan_cpt[:,i])
        # kutta condition
        A[i,1] -= te
        A[i,npan] += te
    end

    # Set RHS
    for i = 1:npan
        RHS[i] += -dot(u_vec, pan_cpt[:,i]) # potential due to freestream
    end

    # Solve linear system
    pan_mu = A\RHS

    # potential on outside of surface
    pot_cpt_out = calc_potential_outer(A_pot_out, pan_mu, pan_cpt, u_vec)
    # potential on inside of surface
    pot_cpt_in = calc_potential_inner(A, pan_mu, pan_cpt, u_vec)
    # freestream potential on boundary collocation points
    pot_freestream = calc_potential_freestream(pan_cpt, u_vec)

    # Approximate velocity at collocation points
    if compute_full_potential
        vel_cpt = calc_vel_from_potential(pot_cpt_out, pan_cpt)
    else
        vel_cpt = calc_vel_from_doublet_strength(pan_mu, pan_cpt)
    end

    # pressure at panels
    pan_pres = calc_pan_pressure(vel_cpt, rho, p0)
    pan_cp = calc_pan_cp(pan_pres, rho, U)

    # force at panels
    pan_force = calc_pan_force(pan_pres, pan_len, pan_norm)

    lift = sum(pan_force[2,:])
    cl = lift/0.5/rho/U^2/chord

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, pan_mu=pan_mu, pan_cpt=pan_cpt, pot_cpt_out=pot_cpt_out, pot_cpt_in=pot_cpt_in, pot_freestream=pot_freestream)
    return result
end

function airfoil_doublet_neumann(pan_vert::Matrix{<:Real}, aoa::Real)
    npan = size(pan_vert,2) - 1

    # Compute panel properties
    pan_cpt = calc_panel_2d_collocation_points(pan_vert) # collocation points
    pan_norm = calc_panel_2d_normals(pan_vert) # panel normals
    pan_tan = calc_panel_2d_tangents(pan_vert) # panel tangents
    pan_len = calc_panel_2d_lengths(pan_vert)

    # Freestream
    U = 1
    alf = deg2rad(aoa)
    chord = 1
    rho = 1
    pstatic = 0
    p0 = pstatic + 0.5*rho*U^2 # stagnation pressure

    u_vec = zeros(2)
    u_vec[1] = cos(alf)
    u_vec[2] = sin(alf)

    A = zeros(npan+1, npan+1) # normal velocity influence coefficient
    B = zeros(npan+1, npan+1) # tangent velocity influence coefficient
    RHS = zeros(npan+1, npan+1)

    # Compute influence coefficients
    for i = 1:npan
        for j = 1:npan
            if i == j
                vel = vel_line_doublet_2d_self(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i], false) # velocity on outside, approach from panel local bottom
            else
                vel = vel_line_doublet_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
            end
            A[i,j] = dot(vel, pan_norm[:,i])
            B[i,j] = dot(vel, pan_tan[:,i])
        end
        # wake panel
        vel = vel_line_doublet_2d(pan_vert[:,1], [1e3;0], pan_cpt[:,i])
        A[i,npan+1] = dot(vel, pan_norm[:,i])
        B[i,npan+1] = dot(vel, pan_tan[:,i])
    end

    # Compute RHS
    for i = 1:npan
        RHS[i] = -dot(u_vec, pan_norm[:,i])
    end

    # Kutta condition
    A[npan+1,1] = 1
    A[npan+1,npan] = -1
    A[npan+1,npan+1] = 1
    RHS[npan+1] = 0

    # remove one panel
    remIdx = Int(round(npan/2))
    keep_idx1 = vcat(collect(1:remIdx-1), collect(remIdx+1:npan+1))
    keep_idx2 = vcat(collect(1:remIdx-1), collect(remIdx+1:npan))
    npan = npan -1

    RHS = RHS[keep_idx1]
    A = A[keep_idx1, keep_idx1]
    B = B[keep_idx1, keep_idx1]
    pan_cpt = pan_cpt[:, keep_idx2]
    pan_norm = pan_norm[:, keep_idx2]
    pan_tan = pan_tan[:, keep_idx2]
    pan_len = pan_len[keep_idx2]

    pan_gam = A\RHS

    # Velocity at collocation points
    vel_cpt = zeros(npan)
    for i = 1:npan
        vel_cpt[i] = dot(B[i,:], pan_gam) + dot(u_vec, pan_tan[:,i])
        # add contribution of eq 11.38 pg 288, 0.5*dmu/dl. Tangent velocity on itself by
        # reconstructing change in potential along airfoil surface despite having constant strength elements
        if i == 1
            l = dist2D(pan_cpt[:,i], pan_cpt[:,i+1])
            vel_cpt[i] += 0.5*(pan_gam[i+1] - pan_gam[i])/l
        elseif i == npan
            l = dist2D(pan_cpt[:,i-1], pan_cpt[:,i])
            vel_cpt[i] += 0.5*(pan_gam[i] - pan_gam[i-1])/l
        else
            l = dist2D(pan_cpt[:,i-1], pan_cpt[:,i+1])
            vel_cpt[i] += 0.5*(pan_gam[i+1] - pan_gam[i-1])/l
        end
    end

    # pressure at panels
    pan_pres = calc_pan_pressure(vel_cpt, rho, p0)
    pan_cp = calc_pan_cp(pan_pres, rho, U)

    # force at panels
    pan_force = calc_pan_force(pan_pres, pan_len, pan_norm)

    lift = sum(pan_force[2,:])
    cl = lift/0.5/rho/U^2/chord

    lift = sum(pan_force[2,:])
    cl = lift/0.5/rho/U^2/chord

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, pan_gam=pan_gam, pan_cpt=pan_cpt)
    return result
end