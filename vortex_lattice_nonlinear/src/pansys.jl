@kwdef struct PanSys{T<:Real, U<:Integer}
    # Panel Mesh
    pan_vert::Matrix{T}
    pan_con::Matrix{U}

    # Panel Properties
    pan_cpt::Matrix{T}
    pan_norm::Matrix{T}
    pan_numpt::Vector{U}
    pan_edgevec::Array{T,3}
    pan_edgelen::Matrix{T}
    pan_edgeuvec::Array{T,3}
    pan_area::Vector{T}
    pan_tang::Array{T,3}
    pan_sin_ti::Matrix{T}
    pan_cos_ti::Matrix{T}
    npan::U
    nvert::U
    pan_gam::Vector{T}
    
    # Panel Neighbors
    pan_neigh_idx::Matrix{U}
    pan_neigh_side::Matrix{U}
    pan_neigh_dir::Matrix{U}
    pan_neigh_num::Vector{U}

    # Airfoil Lookups
    airfoil_idx::Vector{U}
    airfoil_cl::Vector{<:AbstractInterpolation}
    airfoil_cd::Vector{<:AbstractInterpolation}
    airfoil_strips::Matrix{U}
end

function create_pansys(pansys)
    pan_con = pansys.pan_con
    pan_vert = pansys.pan_vert
    npan = size(pan_con, 2)
    nvert = size(pan_vert, 2)
    pan_cpt, pan_norm, pan_numpt, pan_edgevec, pan_edgelen, pan_edgeuvec, pan_area, pan_tang, pan_sin_ti, pan_cos_ti = calc_panprop(pan_vert, pan_con)
    pan_neigh_idx, pan_neigh_side, pan_neigh_dir, pan_neigh_num = calc_neighbors(pan_con, pan_numpt)

    return PanSys(
    # Panel Mesh
    pan_vert = pan_vert,
    pan_con = pan_con,

    # Panel Properties
    pan_cpt = pan_cpt,
    pan_norm = pan_norm,
    pan_numpt = pan_numpt,
    pan_edgevec = pan_edgevec,
    pan_edgelen = pan_edgelen,
    pan_edgeuvec = pan_edgeuvec,
    pan_area = pan_area,
    pan_tang = pan_tang,
    pan_sin_ti = pan_sin_ti,
    pan_cos_ti = pan_cos_ti,
    npan = npan,
    nvert = nvert,
    pan_gam = zeros(npan),
    
    # Panel Neighbors
    pan_neigh_idx = pan_neigh_idx,
    pan_neigh_side = pan_neigh_side,
    pan_neigh_dir = pan_neigh_dir,
    pan_neigh_num = pan_neigh_num,

    # Airfoil Lookups
    airfoil_idx = pansys.airfoil_idx,
    airfoil_cl = pansys.airfoil_cl,
    airfoil_cd = pansys.airfoil_cd,
    airfoil_strips = pansys.airfoil_strips
    )
end