function dist2D(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real})
    # compute 2d distance between two points
    dx = p2[1] - p1[1]
    dz = p2[2] - p1[2]
    return sqrt(dx^2 + dz^2)
end

function dist2Dsq(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real})
    # compute 2d distance squared between two points
    dx = p2[1] - p1[1]
    dz = p2[2] - p1[2]
    return dx^2 + dz^2
end

function norm2D(vec::AbstractVector{<:Real})
    return sqrt(vec[1]^2 + vec[2]^2)
end

function calc_line_tan_2d(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real})
    # compute tangent vector for line element local coordinates
    dx = p2[1] - p1[1]
    dz = p2[2] - p1[2]
    len = sqrt(dx^2 + dz^2)
    return [dx/len; dz/len]
end

function calc_line_norm_2d(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real})
    # compute normal vector for line element local coordinates
    dx = p2[1] - p1[1]
    dz = p2[2] - p1[2]
    len = sqrt(dx^2 + dz^2)
    return [-dz/len; dx/len]
end

function quadrature_transform(a::Real, b::Real, order::Integer)
    # compute scaling factor and transformed integration points
    # for gause quadrature
    gauss_points, weights = gausslobatto(order)
    scale = (b-a)/2
    points = scale.*gauss_points .+ (a+b)/2
    return points, weights, scale
end