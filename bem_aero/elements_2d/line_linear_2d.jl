function pot_line_doublet_linear_2d(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real})
    # Potential induced by linear strength doublet line
    # Potential is expressed in local element frame
    # First output is influence of unit strength at second node
    # Second output is influence of unit strength at second node

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

function pot_line_doublet_linear_2d_self(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real}, limit_plus::Bool=true)
    # Potential induced by linear strength doublet line on itself approached from above
    # Potential is expressed in local element frame
    # First output is influence of unit strength at second node
    # Second output is influence of unit strength at second node

    # Express panel points in local panel frame
    x1 = 0
    x2 = dist2D(p1, p2)
    # x is field point projected onto the local unit vector
    unit_vec = calc_line_tan_2d(p1, p2)
    x = dot(p.-p1, unit_vec)

    phi_a = 0.5*(x/x2-1)
    phi_b = -0.5*(x/x2)

    if !limit_plus
        return -phi_a, -phi_b
    else
        return phi_a, phi_b
    end
end

function pot_line_doublet_linear_2d_int(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real})
    # Potential induced by linear strength doublet line by numerical integration
    # Potential is expressed in local element frame
    # First output is influence of unit strength at second node
    # Second output is influence of unit strength at second node

    n = 5
    vec = p2.-p1 # line between the points
    len = dist2D(p1,p2)
    # ts is coordinate between 0-1 along vec
    ts, weights, scale = quadrature_transform(0, 1, n)
    phia = zero(eltype(p))
    phib = zero(eltype(p))
    for i = 1:n
        p0 = p1 .+ vec.*ts[i]
        pot_point = pot_point_doublet_2d(p0, p, vec)

        gamb = ts[i]
        phib += pot_point*weights[i]*gamb

        gama = (1-ts[i])
        phia += pot_point*weights[i]*gama
    end
    return phia*scale*len, phib*scale*len
end

function pot_line_vortex_linear_2d(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real})
    # Potential induced by linear strength vortex line
    # Potential is expressed in local element frame
    # First output is influence of unit strength at second node
    # Second output is influence of unit strength at second node

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
    r1_2 = r1^2
    r2_2 = r2^2
    theta1 = atan(z,x)
    theta2 = atan(z, x-x2)

    # superposition unit constant and unit ramp
    term1 = -1/(2*pi)*((x-x1)*theta1 - (x-x2)*theta2 + z/2*log(r1_2/r2_2)) # constant part
    term2 = -1/(2*pi)*(x*z/2*log(r1^2/r2^2) + z/2*(x1-x2) + (x^2-x1^2-z^2)/2*theta1 - (x^2-x2^2-z^2)/2*theta2) # ramp part

    phia = term1 - 1/(x2-x1)*term2
    phib = 1/(x2-x1)*term2

    return phia, phib
end

function pot_line_vortex_linear_2d_self(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real}, limit_plus::Bool=true)
    # Potential induced by linear strength vortex line on itself approached from above
    # Potential is expressed in local element frame
    # First output is influence of unit strength at second node
    # Second output is influence of unit strength at second node

    # Express panel points in local panel frame
    x1 = 0
    x2 = dist2D(p1, p2)
    # x is field point projected onto the local unit vector
    unit_vec = calc_line_tan_2d(p1, p2)
    x = dot(p.-p1, unit_vec)

    # superposition unit constant and unit ramp
    term1 = # constant part

    phia = term1 - 1/(x2-x1)*term2
    phib = 1/(x2-x1)*term2

    if !limit_plus
        return -phi_a, -phi_b
    else
        return phi_a, phi_b
    end
end

function pot_line_vortex_linear_2d_int(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real}, use_line_coords::Bool=false)
    # Potential induced by linear strength vortex line by numerical integration
    # Potential is expressed in local element frame
    # First output is influence of unit strength at second node
    # Second output is influence of unit strength at second node

    # line local tangent and normal vectors
    pan_tan = calc_line_tan_2d(p1,p2)
    pan_norm = calc_line_norm_2d(p1,p2)

    if use_line_coords
        # Compute potential in local line coordinate frame. Note that when coordinate frames are rotated, 
        # the potential field is offset by a constant
        # Express field point in local panel frame
        x = dot(p.-p1, pan_tan)
        z = dot(p.-p1, pan_norm)

        # Express panel points in local panel frame
        x1 = 0
        x2 = dist2D(p1, p2)

        p1 = [x1;0]
        p2 = [x2;0]
        p = [x;z]
    end

    n = 5
    vec = p2.-p1 # line between the points
    len = dist2D(p1,p2)
    # ts is coordinate between 0-1 along vec
    ts, weights, scale = quadrature_transform(0, 1, n)
    phia = zero(eltype(p))
    phib = zero(eltype(p))
    for i = 1:n
        p0 = p1 .+ vec.*ts[i]
        pot_point = pot_point_vortex_2d(p0, p)

        gamb = ts[i]
        phib += pot_point*weights[i]*gamb

        gama = (1-ts[i])
        phia += pot_point*weights[i]*gama
    end
    return phia*scale*len, phib*scale*len
end

function vel_line_vortex_linear_2d(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real})
    # Velocity induced by linear strength vortex line
    # Velocity expressed in global frame
    # First output is influence of unit strength at second node
    # Second output is influence of unit strength at second node

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
    # Velocity expressed in global frame
    # First output is influence of unit strength at second node
    # Second output is influence of unit strength at second node

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
    # Velocity expressed in global frame
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
        vel_point = vel_point_vortex_2d(p0, p)

        gamb = ts[i]
        velb .+= vel_point.*weights[i]*(gamb)

        gama = 1-ts[i]
        vela .+= vel_point.*weights[i]*(gama)
    end
    return vela.*scale.*len, velb.*scale.*len
end