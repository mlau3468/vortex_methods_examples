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
    # quit()
    # display(A_pot_out)
    # display(pot_cpt_out)
    # display(pot_cpt_in)

    # Velocity along the surface of airfoil is differentiation of total potential
    vel_cpt = calc_vel_from_potential(pot_cpt_out, pan_cpt)
    # println(sum(vert_gam))
    # display(vert_gam)
    # display(vel_cpt)
    # quit()

    # pressure at panels
    pan_pres = calc_pan_pressure(vel_cpt, rho, p0)
    pan_cp = calc_pan_cp(pan_pres, rho, U)

    # force at panels
    pan_force = calc_pan_force(pan_pres, pan_len, pan_norm)

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
    pot_cpt_out = zeros(npan)
    for i = 1:npan
        pot_cpt_out[i] = dot(A_pot_out[i,:], vert_mu) + dot(u_vec, pan_cpt[:,i])
    end

    # Velocity along the surface of airfoil is differentiation of total potential
    vel_cpt = calc_vel_from_potential(pot_cpt_out, pan_cpt)

    # pressure at panels
    pan_pres = calc_pan_pressure(vel_cpt, rho, p0)
    pan_cp = calc_pan_cp(pan_pres, rho, U)

    # force at panels
    pan_force = calc_pan_force(pan_pres, pan_len, pan_norm)

    lift = sum(pan_force[2,:])
    cl = lift/0.5/rho/U^2/chord

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, vert_mu=vert_mu, pan_cpt=pan_cpt)
    return result
end

function airfoil_vortex_linear_neumann_pot(pan_vert::Matrix{<:Real}, aoa::Real; num_integrate::Bool=false)
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

    # freestream potential on collocation points
    pot_freestream = calc_potential_freestream(pan_cpt, u_vec)

    # derivative of potential is velocity
    #vel_cpt = calc_vel_from_potential(pot_cpt_out, pan_cpt)
    vel_cpt = zeros(npan)
    for i = 1:npan
        vel_cpt[i] = (vert_gam_res[i+1] - vert_gam_res[i])/pan_len[i]
    end

    # display(pot_cpt_in)
    # display(pot_cpt_out)
    display(pot_cpt_in.-pot_cpt_out)
    quit()
end

function airfoil_doublet_neumann_pot(pan_vert::Matrix{<:Real}, aoa::Real)
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

    pan_mu = A\RHS

    # Velocity at collocation points
    vel_cpt = zeros(npan)
    for i = 1:npan
        vel_cpt[i] = dot(B[i,:], pan_mu) + dot(u_vec, pan_tan[:,i])
        # add contribution of eq 11.38 pg 288, 0.5*dmu/dl. Tangent velocity on itself by
        # reconstructing change in potential along airfoil surface despite having constant strength elements
        if i == 1
            l = dist2D(pan_cpt[:,i], pan_cpt[:,i+1])
            vel_cpt[i] += 0.5*(pan_mu[i+1] - pan_mu[i])/l
        elseif i == npan
            l = dist2D(pan_cpt[:,i-1], pan_cpt[:,i])
            vel_cpt[i] += 0.5*(pan_mu[i] - pan_mu[i-1])/l
        else
            l = dist2D(pan_cpt[:,i-1], pan_cpt[:,i+1])
            vel_cpt[i] += 0.5*(pan_mu[i+1] - pan_mu[i-1])/l
        end
    end

    pot_cpt_out = zeros(npan)
    for i = 1:npan
        for j=1:npan
            if i==j
                phi = pot_line_doublet_2d_self(false)*pan_mu[j]
            else
                phi = pot_line_doublet_2d(pan_vert[:,i], pan_vert[:,i+1], pan_cpt[:,i])*pan_mu[j]
            end
            pot_cpt_out[i] += phi
        end
        phi = pot_line_doublet_2d(pan_vert[:,1], [1e3;0], pan_cpt[:,i])*pan_mu[npan+1]
        pot_cpt_out[i] += phi
        pot_cpt_out[i] += dot(u_vec, pan_cpt[:,i])
    end
    vel_cpt = calc_vel_from_potential(pot_cpt_out, pan_cpt)
    display(pot_cpt_out)
    display(vel_cpt)
    quit()

    # pressure at panels
    pan_pres = calc_pan_pressure(vel_cpt, rho, p0)
    pan_cp = calc_pan_cp(pan_pres, rho, U)

    # force at panels
    pan_force = calc_pan_force(pan_pres, pan_len, pan_norm)

    lift = sum(pan_force[2,:])
    cl = lift/0.5/rho/U^2/chord

    lift = sum(pan_force[2,:])
    cl = lift/0.5/rho/U^2/chord

    result = (cl=cl, pan_cp=pan_cp, pan_pres=pan_pres, vel_cpt=vel_cpt, pan_gam=pan_mu, pan_cpt=pan_cpt)
    return result
end