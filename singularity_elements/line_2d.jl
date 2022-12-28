function pot_line_source_2d(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real})
    # pg 234 eq 10.19

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

function pot_line_source_2d_int(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real})
    n = 5
    vec = p2.-p1 # line between the points
    len = dist2D(p1,p2)
    # ts is coordinate between 0-1 along vec
    ts, weights, scale = quadrature_transform(0, 1, n)
    val = 0.0
    for i = 1:n
        p0 = p1 .+ vec.*ts[i]
        val += pot_point_source_2d(p0, p)*weights[i]
    end
    return val*scale*len
end

function vel_line_source_2d(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real}, p::AbstractVector{<:Real})
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

    r1 = dist2D(p1, p)
    r2 = dist2D(p2, p)
    theta1 = atan(z,x)
    theta2 = atan(z, x-x2)

    # Velocity in panel local frame
    upanel = 1/(4*pi)*log(r1^2/r2^2)
    wpanel = 1/(2*pi)*(theta2-theta1)

    # Velocity in global frame
    vel_global = upanel.*pan_tan .+ wpanel.*pan_norm
    return vel_global
end