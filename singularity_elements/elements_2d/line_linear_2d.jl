function vel_line_vortex_linear_2d(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real})
    # Velocity induced by linear strength vortex line

    # line local tangent and normal vectors
    pan_tan = calc_line_tan_2d(p1,p2)
    pan_norm = calc_line_norm_2d(p1,p2)
    
    # Express field point in local panel frame
    x = dot(p.-p1, pan_tan)
    z = dot(p.-p1, pan_norm)

    # Express panel points in local panel frame
    x1 = 0
    x2 = dist2D(p1, p2)

    r1 = dist2D(p1, p)
    r2 = dist2D(p2, p)
    theta1 = atan(z,x)
    theta2 = atan(z, x-x2)

    # Velocity in panel local frame
    upanel_a = z/2/pi*(-1)/(x2-x1)*log(r2/r1) + ((x2-x1)-(x-x1)) /2/pi/(x2-x1) * (theta2-theta1)
    wpanel_a = -((x2-x1)-(x-x1)) /2/pi/(x2-x1) * log(r1/r2) + z/2/pi*(-1)/(x2-x1) * ((x2-x1)/z-(theta2-theta1))
    upanel_b = z/2/pi/(x2-x1)*log(r2/r1) + (x-x1)/2/pi/(x2-x1) * (theta2-theta1)
    wpanel_b = -((x-x1))/2/pi/(x2-x1)*log(r1/r2) + z/2/pi/(x2-x1) * ((x2-x1)/z-(theta2-theta1))

    # Velocity in global frame
    vel_global_a = upanel_a.*pan_tan .+ wpanel_a.*pan_norm
    vel_global_b = upanel_b.*pan_tan .+ wpanel_b.*pan_norm
    
    return vel_global_a, vel_global_b
end

function vel_line_vortex_linear_midpoint_2d(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, limit_plus::Bool=true)
    # Velocity induced by linear strength vortex line at own midpoint approached from above
    # line local tangent and normal vectors
    pan_tan = calc_line_tan_2d(p1,p2)
    pan_norm = calc_line_norm_2d(p1,p2)
    pan_len = dist2D(p1, p2)

    x2 = pan_len
    x = pan_len/2

    # Velocity in panel local frame
    upanel_a = -0.5*(x-x2)/(x2)
    wpanel_a = -1/2/pi
    upanel_b = 0.5*(x)/(x2)
    wpanel_b = 1/2/pi

    if !limit_plus # approach from bottom instead
        upanel_a = -upanel_a
        upanel_b = -upanel_b
    end

    # Velocity in global frame
    vel_global_a = upanel_a.*pan_tan .+ wpanel_a.*pan_norm
    vel_global_b = upanel_b.*pan_tan .+ wpanel_b.*pan_norm
    
    return vel_global_a, vel_global_b
end