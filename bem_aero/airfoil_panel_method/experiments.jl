function airfoil_vortex_linear_dirichlet(pan_vert::Matrix{<:Real}, aoa::Real; num_integrate::Bool=false)
    # set num_integrate = true to use evaluate influence coefficients using numerical quadrature method
    nvert = size(pan_vert, 2)
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

    res1 = airfoil_vortex_linear_neumann(pan_vert, aoa)
    vert_gam_res = res1.vert_gam

    # Compute full potential at each collocation point
    pot_cpt_out = zeros(npan) # potential on outside of surface
    for i = 1:npan
        for j = 1:npan
            if i == j
                phia, phib = pot_line_vortex_linear_2d_self(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i], false)
            else
                phia, phib = pot_line_vortex_linear_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
            end
            pot_cpt_out[i] += vert_gam_res[j]*phia + vert_gam_res[j+1]*phib
        end
        pot_cpt_out[i] += dot(u_vec, pan_cpt[:,i])
    end

    pot_cpt_in = zeros(npan) # potential on inside of surface
    for i = 1:npan
        for j = 1:npan
            if i == j
                phia, phib = pot_line_vortex_linear_2d_self(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
            else
                phia, phib = pot_line_vortex_linear_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
            end
            pot_cpt_in[i] += vert_gam_res[j]*phia + vert_gam_res[j+1]*phib
        end
        pot_cpt_in[i] += dot(u_vec, pan_cpt[:,i])
    end

    pot_freestream = zeros(npan) # freestream potential on boundary collocation points
    for i = 1:npan
        pot_freestream[i] = dot(u_vec, pan_cpt[:,i])
    end
    display(vert_gam_res)
    display(pot_cpt_out)
    # quit()

    # Allocate influence matrices
    A = zeros(nvert, nvert)
    A_pot_out = zeros(npan, nvert)
    RHS = zeros(nvert)

    # compute vortex panel influence coefficients
    for i = 1:npan
        for j = 1:npan
            if i == j 
                phia, phib = pot_line_vortex_linear_2d_self(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
                phia_out, phib_out = pot_line_vortex_linear_2d_self(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i], false)
                # phia = -phia
                # phib = -phib
                # phia_out = -phia_out
                # phib_out = -phib_out
            else
                phia, phib = pot_line_vortex_linear_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
                phia_out = phia
                phib_out = phib
            end

            A[i,j] += phia
            A[i,j+1] += phib
            A_pot_out[i,j] += phia_out
            A_pot_out[i,j+1] += phib_out
        end
    end

    # Set RHS
    for i = 1:npan
        RHS[i] = -dot(u_vec, pan_cpt[:,i]) # potential due to freestream
    end

    # Kutta condition
    A[end, 1] = 1
    A[end, end] = 1
    RHS[end] = 0

    vert_gam = A\RHS
    
    # Compute full potential at each collocation point
    pot_cpt_out = zeros(npan)
    for i = 1:npan
        pot_cpt_out[i] = dot(A_pot_out[i,:], vert_gam) + dot(u_vec, pan_cpt[:,i])
    end
    pot_cpt_in = zeros(npan)
    for i = 1:npan
        pot_cpt_in[i] = dot(A[i,:], vert_gam) + dot(u_vec, pan_cpt[:,i])
    end
    # display(pot_cpt_out)
    # quit()
    pot_cpt = pot_cpt_out

    # Velocity along the surface of airfoil is differentiation of total potential
    vel_vec = zeros(npan-1)
    for i = 1:npan-1
        # finite difference in panel tangent direction
        l = dist2D(pan_cpt[:,i], pan_cpt[:,i+1])
        vel_vec[i] = (pot_cpt_out[i+1] - pot_cpt_out[i])/l
    end

    display(vel_vec)
    quit()

    #Velocity at collocation points
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
    pan_force = zeros(2,npan)
    for i = 1:npan
        pan_force[:,i] .= -pan_pres[i]*pan_len[i] .* pan_norm[:,i]
    end

    lift = sum(pan_force[2,:])
    cl = lift/0.5/rho/U^2/chord

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, vert_gam=vert_gam, pan_cpt=pan_cpt)
    return result
end

function airfoil_vortex_linear_dirichlet2(pan_vert::Matrix{<:Real}, aoa::Real)
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
    A_pot_out = zeros(nvert+1, nvert+1) # Influence coefficeint for potential just outside airfoil
    RHS = zeros(nvert+1)

    pt_add = [0.9;0]
    # Influence coefficients
    for i = 1:npan+1
        if i == npan+1
            cpt = pt_add
        else
            cpt = pan_cpt[:,i]
        end

        for j = 1:npan
            if i == j
                phi_a, phi_b = pot_line_vortex_linear_2d_self(pan_vert[:,j], pan_vert[:,j+1], cpt)
                phi_a_out, phi_b_out = pot_line_vortex_linear_2d_self(pan_vert[:,j], pan_vert[:,j+1], cpt, false)
            else
                phi_a, phi_b = pot_line_vortex_linear_2d(pan_vert[:,j], pan_vert[:,j+1], cpt)
                phi_a_out = phi_a
                phi_b_out = phi_b
            end
            A[i,j] += phi_a
            A[i,j+1] += phi_b

            A_pot_out[i,j] += phi_a_out
            A_pot_out[i,j+1] += phi_b_out
        end

        # wake panel
        A[i,nvert+1] += pot_line_vortex_2d(pan_vert[:,1], [1e3;0], cpt)
        A_pot_out[i,nvert+1] = A[i,nvert+1]
    end

    # Set RHS
    for i = 1:npan
        RHS[i] = -dot(u_vec, pan_cpt[:,i]) # potential due to freestream
    end

    RHS[nvert] = -dot(u_vec, pt_add)

    # Kutta condition
    # A[nvert,1] = 1
    # A[nvert,nvert] = -1
    # A[nvert,nvert+1] = 1

    # A[nvert+1,1] = -1
    # A[nvert+1,2] = 1
    # A[nvert+1,nvert-1] = 1
    # A[nvert+1,nvert] = -1

    # RHS[nvert] = 0
    # RHS[nvert+1] = 0

    # Kutta condition 2
    A[nvert+1,1] = 1
    A[nvert+1,nvert] = -1
    A[nvert+1,nvert+1] = 1
    RHS[nvert+1] = 0

    # Solve linear system
    vert_mu = A\RHS

    # Compute full potential at each collocation point
    pot_cpt = zeros(npan)
    for i = 1:npan
        pot_cpt[i] = dot(A_pot_out[i,:], vert_mu) + dot(u_vec, pan_cpt[:,i])
    end

    # Velocity along the surface of airfoil is differentiation of total potential
    vel_vec = zeros(npan-1)
    for i = 1:npan-1
        # finite difference in panel tangent direction
        l = dist2D(pan_cpt[:,i], pan_cpt[:,i+1])
        vel_vec[i] = (pot_cpt[i+1] - pot_cpt[i])/l
    end

    #Velocity at collocation points
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

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, vert_mu=vert_mu, pan_cpt=pan_cpt)
    return result
end