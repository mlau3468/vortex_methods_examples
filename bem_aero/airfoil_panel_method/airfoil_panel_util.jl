# utility functions for airfoil panel methods

function calc_panel_2d_collocation_points(pan_vert::Matrix{<:Real})
    npan = size(pan_vert,2) - 1
    pan_cpt = zeros(2, npan)
    for i = 1:npan
        pan_cpt[:,i] = midpoint2d(pan_vert[:,i], pan_vert[:,i+1])
    end
    return pan_cpt
end

function calc_panel_2d_tangents(pan_vert::Matrix{<:Real})
    # assume panels are counterclockwise, tangent points in ccw direction of indexing
    npan = size(pan_vert,2) - 1
    pan_tan = zeros(2, npan)
    for i = 1:npan
        pan_tan[:,i] = calc_line_tan_2d(pan_vert[:,i], pan_vert[:,i+1])
    end
    return pan_tan
end

function calc_panel_2d_normals(pan_vert::Matrix{<:Real})
    # assume panels are counterclockwise, normal points outwards of airfoil boundary
    npan = size(pan_vert,2) - 1
    pan_norm = zeros(2, npan)
    for i = 1:npan
        pan_norm[:,i] = calc_line_norm_2d(pan_vert[:,i+1], pan_vert[:,i])
    end
    return pan_norm
end

function calc_panel_2d_lengths(pan_vert::Matrix{<:Real})
    npan = size(pan_vert,2) - 1
    pan_len = zeros(npan)
    for i = 1:npan
        pan_len[i] = dist2D(pan_vert[:,i], pan_vert[:,i+1])
    end
    return pan_len
end

function calc_vel_from_potential(pot_cpt_out::Vector{<:Real}, pan_cpt::Matrix{<:Real})
    npan = size(pan_cpt,2)
    vel_vec = zeros(npan-1)
    for i = 1:npan-1
        # finite difference in panel tangent direction
        l = dist2D(pan_cpt[:,i], pan_cpt[:,i+1])
        vel_vec[i] = (pot_cpt_out[i+1] - pot_cpt_out[i])/l
    end
    vel_cpt = rebuild_vel_cpt(vel_vec)
    return vel_cpt
end

function calc_vel_from_doublet_strength(pan_mu::Vector{<:Real}, pan_cpt::Matrix{<:Real})
    npan = length(pan_mu)
    vel_vec = zeros(npan-1)
    for i = 1:npan-1
        # finite difference in panel tangent direction
        l = dist2D(pan_cpt[:,i], pan_cpt[:,i+1])
        vel_vec[i] = (pan_mu[i+1] - pan_mu[i])/l
    end
    vel_cpt = rebuild_vel_cpt(vel_vec)
    return vel_cpt
end

function rebuild_vel_cpt(vel_vec::Vector{<:Real})
    npan = length(vel_vec) + 1
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
    return vel_cpt
end

function calc_pan_pressure(vel_cpt::Vector{<:Real}, rho::Real, p0::Real)
    npan = length(vel_cpt)
    pan_pres = zeros(npan)
    for i = 1:npan
        pan_pres[i] = p0-0.5*rho*vel_cpt[i]^2
    end
    return pan_pres
end

function calc_pan_cp(pan_pres::Vector{<:Real}, rho::Real, Uinf::Real)
    npan = length(pan_pres)
    pan_cp = zeros(npan)
    for i = 1:npan
        pan_cp[i] = pan_pres[i]/(0.5*rho*Uinf^2)
    end
    return pan_cp
end

function calc_pan_force(pan_pres::Vector{<:Real}, pan_len::Vector{<:Real}, pan_norm::Matrix{<:Real})
    npan = length(pan_pres)
    pan_force = zeros(2, npan)
    for i = 1:npan
        pan_force[:,i] .= -pan_pres[i]*pan_len[i] .* pan_norm[:,i]
    end
    return pan_force
end

function calc_potential_inner(A_pot_in::Matrix{<:Real}, pan_strength::Vector{<:Real}, pan_cpt::Matrix{<:Real}, u_vec::Vector{<:Real})
    npan = length(pan_strength)
    pot_cpt_in = zeros(npan) # potential on inside of surface
    for i = 1:npan
        pot_cpt_in[i] = dot(A_pot_in[i,:], pan_strength) + dot(u_vec, pan_cpt[:,i])
    end
    return pot_cpt_in
end

function calc_potential_inner(A_pot_in::Matrix{<:Real}, B_pot_in::Matrix{<:Real}, pan_strength::Vector{<:Real}, pan_source::Vector{<:Real}, pan_cpt::Matrix{<:Real}, u_vec::Vector{<:Real})
    npan = length(pan_strength)
    pot_cpt_in = zeros(npan) # potential on inside of surface
    for i = 1:npan
        pot_cpt_in[i] = dot(A_pot_in[i,:], pan_strength) + dot(B_pot_in[i,:], pan_source) + dot(u_vec, pan_cpt[:,i])
    end
    return pot_cpt_in
end

function calc_potential_outer(A_pot_out::Matrix{<:Real}, pan_strength::Vector{<:Real}, pan_cpt::Matrix{<:Real}, u_vec::Vector{<:Real})
    npan = length(pan_strength)
    pot_cpt_out = zeros(npan) # potential on outside of surface
    for i = 1:npan
        pot_cpt_out[i] = dot(A_pot_out[i,:], pan_strength) + dot(u_vec, pan_cpt[:,i])
    end
    return pot_cpt_out
end

function calc_potential_outer(A_pot_out::Matrix{<:Real}, B_pot_out::Matrix{<:Real}, pan_strength::Vector{<:Real}, pan_source::Vector{<:Real}, pan_cpt::Matrix{<:Real}, u_vec::Vector{<:Real})
    npan = length(pan_strength)
    pot_cpt_out = zeros(npan) # potential on outside of surface
    for i = 1:npan
        pot_cpt_out[i] = dot(A_pot_out[i,:], pan_strength) + dot(B_pot_out[i,:], pan_source) + dot(u_vec, pan_cpt[:,i])
    end
    return pot_cpt_out
end

function calc_potential_freestream(pan_cpt::Matrix{<:Real}, u_vec::Vector{<:Real})
    npan = size(pan_cpt,2)
    pot_freestream = zeros(npan) # freestream potential on boundary collocation points
    for i = 1:npan
        pot_freestream[i] = dot(u_vec, pan_cpt[:,i])
    end
    return pot_freestream
end