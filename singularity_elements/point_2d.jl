function pot_point_source_2d(p0::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Potential induced by unit strength point source
    # pg 230 eq 10.1
    x = p[1]
    z = p[2]
    x0 = p0[1]
    z0 = p0[2]
    r2 = (x-x0)^2 + (z-z0)^2
    return 1/(2*pi)*log(sqrt(r2))
end

function vel_point_source_2d(p0::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Velocity induced by unit strength point source
    # pg 231 eq 10.2, 10.3
    x = p[1]
    z = p[2]
    x0 = p0[1]
    z0 = p0[2]
    r2 = (x-x0)^2 + (z-z0)^2
    u = 1/(2*pi)*(x-x0)/r2
    w = 1/(2*pi)*(z-x0)/r2
    return [u;w]
end

function vel_point_source_2d_ad(p0::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Velocity induced by unit strength point source by differentiation of potential
    return ForwardDiff.gradient((x) -> pot_point_source_2d(p0, x), p)
end

function pot_point_doublet_2d(p0::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Potential induced by unit strength doublet
    # pg 231 eq 10.4
    x = p[1]
    z = p[2]
    x0 = p0[1]
    z0 = p0[2]
    r2 = (x-x0)^2 + (z-z0)^2
    return -1/(2*pi)*(z-z0)/r2
end

function vel_point_doublet_2d(p0::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Velocity induced by unit strength doublet oriented in z direction
    # pg 231 eq 10.5, 10.6
    x = p[1]
    z = p[2]
    x0 = p0[1]
    z0 = p0[2]
    r2 = (x-x0)^2 + (z-z0)^2
    u = 1/pi*(x-x0)*(z-z0)/r2^2
    w = -1/(2*pi)* ((x-x0)^2 - (z-z0)^2) /r2^2
    return [u;w]
end

function vel_point_doublet_2d_ad(p0::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Velocity induced by unit strength doublet oriented in z direction by differentiation of potential
    return ForwardDiff.gradient((x) -> pot_point_doublet_2d(p0, x), p)
end

function pot_point_vortex_2d(p0::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Potential induced by unit strength vortex
    # pg 232 eq. 10.8
    x = p[1]
    z = p[2]
    x0 = p0[1]
    z0 = p0[2]
    return -1/(2*pi)*atan(z-z0, x-x0)
end

function vel_point_vortex_2d(p0::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Velocity induced by unit strength vortex
    # pg 232 eq 10.9, 10.10
    x = p[1]
    z = p[2]
    x0 = p0[1]
    z0 = p0[2]
    r2 = (x-x0)^2 + (z-z0)^2
    u = 1/(2*pi)*(z-z0)/r2
    w = -1/(2*pi)*(x-x0)/r2
    return [u;w]
end

function vel_point_vortex_2d_ad(p0::AbstractVector{T}, p::AbstractVector{T}) where {T<:Real}
    # Velocity induced by unit strength vortex by differentiation of potential
    return ForwardDiff.gradient((x) -> pot_point_vortex_2d(p0, x), p)
end