function pot_line_doublet_linear_2d(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real})
    # Potential induced by linear strength doublet line

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

    # phi_a = -1/(2*pi) * (theta2 - theta1 - 1/(x2-x1) * (x*(theta2-theta1) + z/2*log(r2^2/r1^2)))
    # phi_b = -1/(2*pi*(x2-x1)) * (x*(theta2-theta1) + z/2*log(r2^2/r1^2))
    phi_a = 1/(2*pi) * (x/x2*(theta2-theta1) + z/x2*log(r2/r1) - (theta2-theta1))
    phi_b = -1/(2*pi) * (x/x2*(theta2-theta1) + z/x2*log(r2/r1))
    return phi_a, phi_b
end

function pot_line_doublet_linear_2d_self(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, limit_plus::Bool=true)
    # Potential induced by linear strength doublet line on itself approached from above

    # Express panel points in local panel frame
    x1 = 0
    x2 = dist2D(p1, p2)
    x = x2/2

    phi_a = 0.5*(x/x2-1)
    phi_b = -0.5*(x/x2)

    if !limit_plus
        return -phi_a, -phi_b
    else
        return phi_a, phi_b
    end
end

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

function vel_line_vortex_linear_2d_int(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real})
    # Velocity induced by unit linear strength vortex line by numerical integration
    # First output is influence of unit trength at second node
    # Second output is influence of unit trength at second node
    n = 5
    vec = p2.-p1 # line between the points
    len = dist2D(p1,p2)
    # ts is coordinate between 0-1 along vec
    ts, weights, scale = quadrature_transform(0, 1, n)
    vela = zeros(eltype(p),2)
    velb = zeros(eltype(p),2)
    for i = 1:n
        p0 = p1 .+ vec.*ts[i]
        val = vel_point_vortex_2d(p0, p)

        gamb = ts[i]
        velb .+= val.*weights[i]*(gamb)

        gama = 1-ts[i]
        vela .+= val.*weights[i]*(gama)
    end
    return vela.*scale.*len, velb.*scale.*len
end