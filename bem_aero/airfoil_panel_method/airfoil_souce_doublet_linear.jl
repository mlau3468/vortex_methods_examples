
function airfoil_source_doublet_linear_dirichlet(pan_vert::Matrix{<:Real}, aoa::Real)
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

    A = zeros(nvert+1, nvert+1) # Influence coefficents for potential just inside airfoil
    A_pot_out = zeros(nvert+1, nvert+1) # Influence coefficents for potential just outside airfoil
    B = zeros(npan, npan) # source panel influence coefficents
    RHS = zeros(nvert+1)

    # Influence coefficients
    for i = 1:npan
        for j = 1:npan
            if i == j
                phi_a, phi_b = pot_line_doublet_linear_2d_self(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
                phi_a_out, phi_b_out = pot_line_doublet_linear_2d_self(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i], false)
            else
                phi_a, phi_b = pot_line_doublet_linear_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
                phi_a_out = phi_a
                phi_b_out = phi_b
            end
            A[i,j] += phi_a
            A[i,j+1] += phi_b

            A_pot_out[i,j] += phi_a_out
            A_pot_out[i,j+1] += phi_b_out
        end
        # wake panel
        A[i,nvert+1] += pot_line_doublet_2d(pan_vert[:,1], [1e3;0], pan_cpt[:,i])
        A_pot_out[i,nvert+1] = A[i,nvert+1]
    end

    # Source panel influence coefficeints
    for i = 1:npan
        for j = 1:npan
            B[i,j] = pot_line_source_2d(pan_vert[:,j], pan_vert[:,j+1],  pan_cpt[:,i])
        end
    end

    # Source strengths
    # Source strengths selected such that the internal velocity potential is equal to the 
    # freestream velocity potential
    pan_source = zeros(npan)
    for i = 1:npan
        pan_source[i] = dot(calc_line_norm_2d(pan_vert[:,i], pan_vert[:,i+1]), u_vec)
    end

    #Set RHS
    for i = 1:npan
        RHS[i] = -dot(B[i,:], pan_source) # source panel
    end

    # Kutta condition
    A[nvert,1] = 1
    A[nvert,nvert] = -1
    A[nvert,nvert+1] = 1

    A[nvert+1,1] = -1
    A[nvert+1,2] = 1
    A[nvert+1,nvert-1] = 1
    A[nvert+1,nvert] = -1

    RHS[nvert] = 0
    RHS[nvert+1] = 0

    # Solve linear system
    vert_mu = A\RHS

    # potential on outside of surface
    pot_cpt_out = calc_potential_outer(A_pot_out, B, vert_mu, pan_source, pan_cpt, u_vec)
    # potential on inside of surface
    pot_cpt_in = calc_potential_inner(A, B, vert_mu, pan_source, pan_cpt, u_vec)
    # freestream potential on boundary collocation points
    pot_freestream = calc_potential_freestream(pan_cpt, u_vec)

    # Velocity along the surface of airfoil is differentiation of total potential
    vel_cpt = calc_vel_from_potential(pot_cpt_out, pan_cpt)

    # The following results in a noisy derivative
    # vel_cpt2 = zeros(npan)
    # for i = 1:npan
    #     # dphia, dphib = d_pot_line_doublet_linear_2d_self(pan_vert[:,i], pan_vert[:,i+1], pan_cpt[:,i], false)
    #     #vel_cpt2[i] = (dphia*vert_mu[i] + dphib*vert_mu[i+1]) + dot(u_vec, pan_tan[:,i])
    #     vel_cpt2[i] = (vert_mu[i+1] - vert_mu[i])/pan_len[i] + dot(u_vec, pan_tan[:,i])
    # end

    # pressure at panels
    pan_pres = calc_pan_pressure(vel_cpt, rho, p0)
    pan_cp = calc_pan_cp(pan_pres, rho, U)

    # force at panels
    pan_force = calc_pan_force(pan_pres, pan_len, pan_norm)

    lift = sum(pan_force[2,:])
    cl = lift/0.5/rho/U^2/chord

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, vert_mu=vert_mu, pan_cpt=pan_cpt, pot_cpt_out=pot_cpt_out, pot_cpt_in=pot_cpt_in, pot_freestream=pot_freestream)
    return result
end