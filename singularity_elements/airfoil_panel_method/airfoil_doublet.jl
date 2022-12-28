function airfoil_doublet_dirichlet(pan_vert::Matrix{<:Real}, aoa::Real)
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

    A = zeros(npan, npan)
    RHS = zeros(npan)

    # Influence coefficients of doublet
    for i = 1:npan
        for j = 1:npan
            if i == j
                A[i,j] = pot_line_doublet_2d_self() # limit approaching from inside airfoil
            else
                A[i,j] = pot_line_doublet_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
            end
        end
        # trailing edge influence
        te = pot_line_doublet_2d(pan_vert[:,1], [1e3;0.0], pan_cpt[:,i])
        A[i,1] -= te
        A[i,npan] += te
    end

     # Set RHS
     for i = 1:npan
        RHS[i] += -dot(u_vec, pan_cpt[:,i]) # potential due to freestream
    end

    # Solve linear system
    pan_mu = A\RHS

    # Velocity along the surface of airfoil is differentiation of total potential
    vel_vec = zeros(npan-1)
    for i = 1:npan-1
        # finite difference in panel tangent direction
        l = dist2D(pan_cpt[:,i], pan_cpt[:,i+1])
        vel_vec[i] = (pan_mu[i+1] - pan_mu[i])/l
    end

    # Velocity at collocation points
    vel_cpt = zeros(npan)
    for i = 1:npan
        if i == 1
            vel_cpt[i] = vel_vec[i]
        elseif i == npan
            vel_cpt[i] = vel_vec[i-1]
        else
            vel_cpt[i] = (vel_vec[i-1] + vel_vec[i])/2
        end
    end

    # pressure at panels
    pan_pres = zeros(npan)
    pan_cp = zeros(npan)
    for i = 1:npan
        pan_pres[i] = p0-0.5*rho*vel_cpt[i]^2
        pan_cp[i] = pan_pres[i]/(0.5*rho*U^2)
    end

    # force at panels
    pan_force = zeros(2, npan)
    for i = 1:npan
        pan_force[:,i] .= -pan_pres[i]*pan_len[i] .* pan_norm[:,i]
    end

    lift = sum(pan_force[2,:])
    cl = lift/0.5/rho/U^2/chord

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, pan_mu=pan_mu, pan_cpt=pan_cpt)
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
    RHS = zeros(npan+1)

    for i = 1:npan
        for j = 1:npan
            if i == j
                vel = vel_line_doublet_2d_midpoint(pan_vert[:,j+1], pan_vert[:,j])
            else
                vel = -vel_line_doublet_2d(pan_vert[:,j+1], pan_vert[:,j], pan_cpt[:,i])
            end
            A[i,j] = dot(vel, pan_norm[:,i])
            B[i,j] = dot(vel, pan_tan[:,i])
        end
        # trailing edge
        # vel = vel_line_doublet_2d([1e3;0], pan_vert[:,1], pan_cpt[:,i])
        vel = -vel_line_doublet_2d(pan_vert[:,1], [1e3;0] , pan_cpt[:,i])
        A[i,npan+1] = dot(vel, pan_norm[:,i])
        B[i,npan+1] = dot(vel, pan_tan[:,i])
    end

    for i = 1:npan
        RHS[i] = -dot(u_vec, pan_norm[:,i])
    end

    # Kutta condition
    A[npan+1,1] = 1
    A[npan+1,npan] = -1
    A[npan+1,npan+1] = 1
    RHS[npan+1] = 0

    # Solve linear system
    pan_mu = A\RHS

    # Velocity at collocation points
    vel_cpt = zeros(npan)
    for i = 1:npan
        vel_cpt[i] = dot(B[i,:], pan_mu) + dot(u_vec, pan_tan[:,i])
    end

    # pressure at panels
    pan_pres = zeros(npan)
    pan_cp = zeros(npan)
    for i = 1:npan
        pan_pres[i] = p0-0.5*rho*vel_cpt[i]^2
        pan_cp[i] = pan_pres[i]/(0.5*rho*U^2)
    end

    # force at panels
    pan_force = zeros(2, npan)
    for i = 1:npan
        pan_force[:,i] .= -pan_pres[i]*pan_len[i] .* pan_norm[:,i]
    end

    lift = sum(pan_force[2,:])
    cl = lift/0.5/rho/U^2/chord

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, pan_mu=pan_mu, pan_cpt=pan_cpt)
    return result
end