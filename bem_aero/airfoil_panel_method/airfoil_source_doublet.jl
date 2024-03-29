function airfoil_sourcedoublet_dirichlet(pan_vert::Matrix{<:Real}, aoa::Real; compute_full_potential::Bool=false)
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

    A = zeros(npan+1, npan+1) # Influence coefficent for potential just inside airfoil
    A_pot_out = zeros(npan+1, npan+1) # Influence coefficeint for potential just outside airfoil
    B = zeros(npan, npan)
    RHS = zeros(npan+1)

    # Influence coefficients of doublet and sources
    for i = 1:npan
        for j = 1:npan
            if i == j
                A[i,j] = pot_line_doublet_2d_self() # limit approaching from inside airfoil
                A_pot_out[i,j] = pot_line_doublet_2d_self(false) # limit approaching from outside airfoil
            else
                A[i,j] = pot_line_doublet_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
                A_pot_out[i,j] = A[i,j]
            end
            B[i,j] = pot_line_source_2d(pan_vert[:,j], pan_vert[:,j+1],  pan_cpt[:,i])
        end
        # trailing edge influence
        A[i,npan+1] = pot_line_doublet_2d(pan_vert[:,1], [1e3;0.0], pan_cpt[:,i])
    end

    # Source strengths
    # Source strengths selected such that the internal velocity potential is equal to the 
    # freestream velocity potential
    pan_source = zeros(npan)
    for i = 1:npan
        pan_source[i] = dot(calc_line_norm_2d(pan_vert[:,i], pan_vert[:,i+1]), u_vec)
    end

    # Set RHS
    for i = 1:npan
        RHS[i] = -dot(B[i,:], pan_source) # source panel
    end

    # Kutta condition
    A[npan+1,1] = 1
    A[npan+1,npan] = -1
    A[npan+1,npan+1] = 1
    RHS[npan+1] = 0

    # Solve linear system
    pan_mu = A\RHS

    # potential on outside of surface
    pot_cpt_out = calc_potential_outer(A_pot_out, B, pan_mu, pan_source, pan_cpt, u_vec)
    # potential on inside of surface
    pot_cpt_in = calc_potential_inner(A, B, pan_mu, pan_source, pan_cpt, u_vec)
    # freestream potential on boundary collocation points
    pot_freestream = calc_potential_freestream(pan_cpt, u_vec)

    # Velocity along the surface of airfoil is differentiation of potential
    vel_cpt = calc_vel_from_potential(pot_cpt_out, pan_cpt)

    # Alternatively, the velocity can be found by first considering the contribution of the perturbation potential due to the 
    # doublets, computed by the line below. 
    use_below = false
    if use_below
        vel_cpt = calc_vel_from_const_doublet_panels(pan_mu, pan_cpt)
        # As the perturbation was performed on top of the freestream potential,
        # the velocity contirbution due to the freestream potential is added.
        for i = 1:npan
            vel_cpt[i] += dot(u_vec, pan_tan[:,i]) # note that the freestream velocity vectors is by definition the derivative of the freestream potential
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

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, pan_mu=pan_mu, pan_source=pan_source, pan_cpt=pan_cpt, pot_cpt_out=pot_cpt_out, pot_cpt_in=pot_cpt_in, pot_freestream=pot_freestream)
    return result
end