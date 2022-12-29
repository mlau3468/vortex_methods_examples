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

    A = zeros(npan, npan) # Influence coefficent for potential just inside airfoil
    A_pot_out = zeros(npan, npan) # Influence coefficeint for potential just outside airfoil
    B = zeros(npan, npan)
    RHS = zeros(npan)

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
        # wake panel
        te = pot_line_doublet_2d(pan_vert[:,1], [1e3;0.0] , pan_cpt[:,i])
        # kutta condition
        A[i,1] -= te
        A[i,npan] += te
    end

    # Source strengths
    pan_source = zeros(npan)
    for i = 1:npan
        pan_source[i] = dot(calc_line_norm_2d(pan_vert[:,i], pan_vert[:,i+1]), u_vec)
    end

    # Set RHS
    for i = 1:npan
        RHS[i] = -dot(B[i,:], pan_source) # source panel
    end

    # Solve linear system
    pan_mu = A\RHS

    # Compute full potential at each collocation point
    pot_cpt_out = zeros(npan)  # potential on outside of surface
    for i = 1:npan
        pot_cpt_out[i] = dot(A_pot_out[i,:], pan_mu) + dot(B[i,:], pan_source) + dot(u_vec, pan_cpt[:,i])
    end
    pot_cpt_in = zeros(npan) # potential on inside of surface
    for i = 1:npan
        pot_cpt_in[i] = dot(A[i,:], pan_mu) + dot(B[i,:], pan_source) + dot(u_vec, pan_cpt[:,i])
    end
    pot_freestream = zeros(npan) # potential on boundary collocation points
    for i = 1:npan
        pot_freestream[i] = dot(u_vec, pan_cpt[:,i])
    end

    # Velocity along the surface of airfoil is differentiation of potential
    # note dot(u_vec, pan_tan[:,i]) is needed as this is a pertubation potential formulation
    # dot(u_vec, pan_tan[:,i]) includes the effect of the gradient of freestream potential function
    vel_vec = zeros(npan-1)
    for i = 1:npan-1
        # finite difference in panel tangent direction
        l = dist2D(pan_cpt[:,i], pan_cpt[:,i+1])
        if compute_full_potential
            vel_vec[i] = (pot_cpt_out[i+1] - pot_cpt_out[i])/l
        else
            vel_vec[i] = (pan_mu[i+1] - pan_mu[i])/l + dot(u_vec, pan_tan[:,i])
        end
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

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, pan_mu=pan_mu, pan_source=pan_source, pan_cpt=pan_cpt, pot_cpt_out=pot_cpt_out, pot_cpt_in=pot_cpt_in, pot_freestream=pot_freestream)
    return result
end