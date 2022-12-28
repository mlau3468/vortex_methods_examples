function airfoil_vortex_linear_neumann(pan_vert::Matrix{<:Real}, aoa::Real)
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
                uw_a, uw_b = vel_line_vortex_linear_midpoint_2d(pan_vert[:,j], pan_vert[:,j+1], false)
            else
                uw_a, uw_b = vel_line_vortex_linear_2d(pan_vert[:,j], pan_vert[:,j+1], pan_cpt[:,i])
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

    # velocity at panels
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