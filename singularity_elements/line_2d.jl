function pot_line_source_2d(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Potential induced by unit strength source line
    # pg 234 eq 10.19 corrected suing appendix B.11

    # line local tangent and normal vectors
    pan_tan = calc_line_tan_2d(p1,p2)
    pan_norm = calc_line_norm_2d(p1,p2)

    # Express field point in local panel frame
    x = dot(p.-p1, pan_tan)
    z = dot(p.-p1, pan_norm)

    # Express panel points in local panel frame
    x1 = 0
    x2 = dist2D(p1, p2)

    term1 = (x-x1)*log((x-x1)^2 + z^2)
    term2 = (x-x2)*log((x-x2)^2 + z^2)
    term3 = atan(z, x-x2) - atan(z, x-x1)

    phi = 1/(4*pi)*(term1 - term2 -2*(x2-x1) + 2*z*term3)
    return phi
end

function pot_line_source_2d_int(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Potential induced by unit strength source line by numerical integration
    n = 5
    vec = p2.-p1 # line between the points
    len = dist2D(p1,p2)
    # ts is coordinate between 0-1 along vec
    ts, weights, scale = quadrature_transform(0, 1, n)
    val = zero(T)
    for i = 1:n
        p0 = p1 .+ vec.*ts[i]
        val += pot_point_source_2d(p0, p)*weights[i]
    end
    return val*scale*len
end

function vel_line_source_2d(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Velocity induced by unit strength source line
    # pg 234 eq 10.20, 10,21

    # line local tangent and normal vectors
    pan_tan = calc_line_tan_2d(p1,p2)
    pan_norm = calc_line_norm_2d(p1,p2)

    # Express field point in local panel frame
    x = dot(p.-p1, pan_tan)
    z = dot(p.-p1, pan_norm)

    # Express panel points in local panel frame
    x1 = 0
    x2 = dist2D(p1, p2)

    r1_2 = dist2Dsq(p1, p)
    r2_2 = dist2Dsq(p2, p)
    theta1 = atan(z,x)
    theta2 = atan(z, x-x2)

    # Velocity in panel local frame
    upanel = 1/(4*pi)*log(r1_2/r2_2)
    wpanel = 1/(2*pi)*(theta2-theta1)

    # Velocity in global frame
    vel_global = upanel.*pan_tan .+ wpanel.*pan_norm
    return vel_global
end

function vel_line_source_2d_int(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Velocity induced by unit strength source line by numerical integration
    n = 5
    vec = p2.-p1 # line between the points
    len = dist2D(p1,p2)
    # ts is coordinate between 0-1 along vec
    ts, weights, scale = quadrature_transform(0, 1, n)
    val = zeros(T,2)
    for i = 1:n
        p0 = p1 .+ vec.*ts[i]
        val .+= vel_point_source_2d(p0, p)*weights[i]
    end
    return val.*scale.*len
end

function pot_line_doublet_2d(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Potential induced by unit doublet line
    # eq 10.28 pg 235

    # line local tangent and normal vectors
    pan_tan = calc_line_tan_2d(p1,p2)
    pan_norm = calc_line_norm_2d(p1,p2)

    # Express field point in local panel frame
    x = dot(p.-p1, pan_tan)
    z = dot(p.-p1, pan_norm)

    # Express panel points in local panel frame
    x1 = 0
    x2 = dist2D(p1, p2)

    theta1 = atan(z,x)
    theta2 = atan(z, x-x2)

    return -1/(2*pi)*(theta2-theta1)
end

function pot_line_doublet_2d_int(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
end

function vel_line_doublet_2d(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
end

function vel_line_doublet_2d_int(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
end

function pot_line_vortex_2d(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Potential induced by unit strength vortex line
    # eq 10.37 pg 237
    # line local tangent and normal vectors
    pan_tan = calc_line_tan_2d(p1,p2)
    pan_norm = calc_line_norm_2d(p1,p2)

    # Express field point in local panel frame
    x = dot(p.-p1, pan_tan)
    z = dot(p.-p1, pan_norm)

    # Express panel points in local panel frame
    x1 = 0
    x2 = dist2D(p1, p2)

    r1_2 = dist2Dsq(p1, p)
    r2_2 = dist2Dsq(p2, p)
    theta1 = atan(z,x)
    theta2 = atan(z, x-x2)

    return -1/(2*pi)*((x-x1)*theta1 - (x-x2)*theta2 + z/2*log(r1_2/r2_2))
end

function pot_line_vortex_2d_int(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Potential induced by unit strength vortex line by numerical integration
    n = 5
    vec = p2.-p1 # line between the points
    len = dist2D(p1,p2)
    # ts is coordinate between 0-1 along vec
    ts, weights, scale = quadrature_transform(0, 1, n)
    val = zero(T)
    for i = 1:n
        p0 = p1 .+ vec.*ts[i]
        val += pot_point_vortex_2d(p0, p)*weights[i]
    end
    return val*scale*len
end

function vel_line_vortex_2d(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Velocity induced by unit strength vortex line
    # line local tangent and normal vectors
    pan_tan = calc_line_tan_2d(p1,p2)
    pan_norm = calc_line_norm_2d(p1,p2)

    # Express field point in local panel frame
    x = dot(p.-p1, pan_tan)
    z = dot(p.-p1, pan_norm)

    # Express panel points in local panel frame
    x1 = 0
    x2 = dist2D(p1, p2)

    r1_2 = dist2Dsq(p1, p)
    r2_2 = dist2Dsq(p2, p)
    theta1 = atan(z,x)
    theta2 = atan(z, x-x2)

    # Velocity in panel local frame
    upanel = 1/(2*pi)*(theta2-theta1)
    wpanel = -1/(4*pi)*log(r1_2/r2_2)

    # Velocity in global frame
    vel_global = upanel.*pan_tan .+ wpanel.*pan_norm
    return vel_global
end

function vel_line_vortex_2d_int(p1::AbstractVector{T}, p2::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Velocity induced by unit strength vortex line by numerical integration
    n = 5
    vec = p2.-p1 # line between the points
    len = dist2D(p1,p2)
    # ts is coordinate between 0-1 along vec
    ts, weights, scale = quadrature_transform(0, 1, n)
    val = zeros(T,2)
    for i = 1:n
        p0 = p1 .+ vec.*ts[i]
        val .+= vel_point_vortex_2d(p0, p)*weights[i]
    end
    return val.*scale.*len
end