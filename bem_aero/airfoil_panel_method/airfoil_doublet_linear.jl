function airfoil_doublet_linear_dirichlet(pan_vert::Matrix{<:Real}, aoa::Real)
    nvert = size(pan_vert,2)
    npan = nvert - 1

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

    A = zeros(nvert+1, nvert+1) # Influence coefficent for potential just inside airfoil
    # A_pot_out = zeros(nvert+1, nvert+1) # Influence coefficeint for potential just outside airfoil
    RHS = zeros(nvert+1)

    # Influence coefficients
    for i = 1:npan
        for j = 1:npan
            if i == j
                phi_a, phi_b = pot_line_doublet_linear_2d_self(pan_vert[:,j], pan_vert[:,j+1], false)
            else
                phi_a, phi_b = pot_line_doublet_linear_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
            end
            A[i,j] += phi_a
            A[i,j+1] += phi_b
        end
        # trailing edge
        A[i,nvert+1] = pot_line_doublet_2d(pan_vert[:,1], [1e3;0], pan_cpt[:,i])
    end

    # Set RHS
    for i = 1:npan
        RHS[i] += -dot(u_vec, pan_cpt[:,i]) # potential due to freestream
    end

    # Kutta condition
    A[nvert,1] = 1
    A[nvert,nvert] = -1
    A[nvert,nvert+1] = 1

    A[nvert+1,1] = -1
    A[nvert+1,2] = 1
    A[nvert+1,nvert-1] = 1
    A[nvert+1,nvert] = -1
    RHS[nvert+1] = 0

    display(A)
    
    # Solve linear system
    vert_mu = A\RHS

    # Velocity at collocation points
    vel_cpt = zeros(npan)
    for i = 1:npan
        vel_cpt[i] = (vert_mu[i+1] - vert_mu[i])/pan_len[i]
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

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, vert_mu=vert_mu, pan_cpt=pan_cpt)
    return result
end