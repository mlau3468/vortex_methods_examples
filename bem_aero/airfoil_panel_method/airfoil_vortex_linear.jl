function airfoil_vortex_linear_neumann(pan_vert::Matrix{<:Real}, aoa::Real; num_integrate::Bool=false)
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

    # Allocate influence matrices
    A = zeros(nvert, nvert) # normal velocity influence coefficient
    B = zeros(npan, nvert) # tangent velocity influence coefficient
    RHS = zeros(nvert)

    # compute vortex panel influence coefficients
    for i = 1:npan
        for j = 1:npan
            if i == j 
                uw_a, uw_b = vel_line_vortex_linear_2d_self(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i], false)
            else
                if num_integrate
                    uw_a, uw_b = vel_line_vortex_linear_2d_int(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
                else
                    uw_a, uw_b = vel_line_vortex_linear_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
                end
            end
            A[i,j] += dot(uw_a, pan_norm[:,i])
            A[i,j+1] += dot(uw_b, pan_norm[:,i])
            B[i,j] += dot(uw_a, pan_tan[:,i])
            B[i,j+1] += dot(uw_b, pan_tan[:,i])
        end
    end

    # compute RHS
    for i = 1:npan
        # Freestream contribution to RHS 
        RHS[i] = -dot(u_vec, pan_norm[:,i])
    end

    # Kutta condition
    A[end, 1] = 1
    A[end, end] = 1
    RHS[end] = 0

    vert_gam = A\RHS

    # velocity at collocation points
    vel_cpt = zeros(npan)
    for i = 1:npan
        vel_cpt[i] = dot(B[i,:], vert_gam) + dot(u_vec, pan_tan[:,i])
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

    # res1 = airfoil_vortex_linear_neumann(pan_vert, aoa)
    # vert_gam_res = res1.vert_gam

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
            else
                phia, phib = pot_line_vortex_linear_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
                phia_out = phia
                phib_out = phib
            end

            A[i,j] += phia
            A[i,j+1] += phib
            A_pot_out[i,j] += phia_out
            A_pot_out[i,j] += phib_out
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
    display(A)
    display(vert_gam)
    quit()

    # Compute full potential at each collocation point
    pot_cpt = zeros(npan)
    for i = 1:npan
        pot_cpt[i] = dot(A_pot_out[i,:], vert_gam) + dot(u_vec, pan_cpt[:,i])
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
    pan_force = zeros(2,npan)
    for i = 1:npan
        pan_force[:,i] .= -pan_pres[i]*pan_len[i] .* pan_norm[:,i]
    end

    lift = sum(pan_force[2,:])
    cl = lift/0.5/rho/U^2/chord

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, vert_gam=vert_gam, pan_cpt=pan_cpt)
    return result
end