"""
    calc_vortlatsoln!(pan_pres, pan_force, pan_gam, pan_dgamdt, pan_neigh_idx, pan_velwake, pan_velcpt, pan_edgeuvec, pan_edgelen, pan_area, pan_norm, npan, rhoinf, uinf, elem_type)
Computes vortex lattice pressure and forces. Mutates pan_pres and pan_force. 
"""
function calc_vortlatsoln!(pansys::PanSys{T,U}, rhoinf::Real, uinf::Vector{<:Real}) where {T<:Real,U<:Integer}
    pan_area = pansys.pan_area
    pan_edgelen = pansys.pan_edgelen
    pan_norm = pansys.pan_norm
    pan_neigh_idx = pansys.pan_neigh_idx
    pan_edgeuvec = pansys.pan_edgeuvec
    pan_gam = pansys.pan_gam
    npan = pansys.npan
    pan_velwake = pansys.pan_velwake
    pan_pres = pansys.pan_pres
    pan_force = pansys.pan_force

    @views @inbounds for i = 1:npan
        val = 0.0
        # direction of increasing column, towards te
        if pan_neigh_idx[4,i] > 0
            gam2 = pan_gam[pan_neigh_idx[4,i]]
        else
            gam2 = 0.0
        end 
        val = val .+ dot(uinf.- pan_velcpt[:,i]+pan_velwake[:,i], pan_edgeuvec[:,1,i]).*(pan_gam[i]-gam2)./pan_edgelen[1,i]
        # direction of increasing rows, perpendicular to te
        if pan_neigh_idx[1,i] > 0
            gam2 = pan_gam[pan_neigh_idx[1,i]]
        else
            gam2 = 0.0
        end
        val = val .+ dot(uinf.- pan_velcpt[:,i]+pan_velwake[:,i], pan_edgeuvec[:,4,i]).*(pan_gam[i]-gam2)./pan_edgelen[4,i]
        pan_pres[i] = -rhoinf*val
        pan_force[:,i] = pan_pres[i].*pan_area[i].*pan_norm[:,i]
    end
    return
end

"""
    vel_line!(v_dou, loc, gam, pts, edgeUni, edgeLen)
Velocity induced by constant 1 intensity vortex line.
Linear core (Rankine vortex) regularization.
# Arguments
- `v_dou::Vector`: vector to store velocity vector in
- `pts::Array`: point cloud of the lifting line system
- `pts_idx::Vector`: The 2 points that define the line
- `edge_uni::Vector`: unit vector of this line
- `edge_len::Real`: length of line
- `loc::Vector`: location to calculate velocity at
- `gam::Real`: circulation intensity
"""
function vel_line!(v_dou, pts, pts_idx, edge_uni, edge_len, loc, gam, r_rankine, r_cutoff)
    #av = loc.-pts[:,1]
    avx = loc[1] - pts[1,pts_idx[1]]
    avy = loc[2] - pts[2,pts_idx[1]]
    avz = loc[3] - pts[3,pts_idx[1]]

    #ai = dot(av,edge_uni)
    ai = avx*edge_uni[1] + avy*edge_uni[2] + avz*edge_uni[3]

    #R1 = norm(av)
    R1 = sqrt(avx^2 + avy^2 + avz^2)

    #R2 = norm(loc.-pts[:,pts_idx[2]])
    R2 = sqrt((loc[1]-pts[1,pts_idx[2]])^2 + (loc[2]-pts[2,pts_idx[2]])^2 + (loc[3]-pts[3,pts_idx[2]])^2)

    #hv = av .- ai.*edgeUni
    hvx = avx - ai*edge_uni[1]
    hvy = avy - ai*edge_uni[2]
    hvz = avz - ai*edge_uni[3]

    #hi = norm(hv)
    hi = sqrt(hvx^2 + hvy^2 + hvz^2)

    if hi > edge_len.*r_rankine
        #vdou = ((edge_len.-ai)./R2 .+ai./R1) ./ (hi.^2) .* cross(edge_uni, hv)
        temp = ((edge_len.-ai)./R2 .+ ai./R1) ./ (hi.^2)
        v_dou[1] += temp * (edge_uni[2]*hvz-edge_uni[3]*hvy)* gam/(4*pi)
        v_dou[2] += temp * (-edge_uni[1]*hvz+edge_uni[3]*hvx)* gam/(4*pi)
        v_dou[3] += temp * (edge_uni[1]*hvy-edge_uni[2]*hvx)* gam/(4*pi)
    else
        if R1 > edge_len.*r_cutoff && R2 > edge_len.*r_cutoff
            r_ran = r_rankine .* edge_len
            #vdou = ((edgeLen.-ai)./R2 .+ai./R1) ./ (r_ran.^2) .* cross(edgeUni, hv)
            temp = ((edge_len.-ai)./R2 + ai./R1) ./ (r_ran.^2)
            v_dou[1] += temp * (edge_uni[2]*hvz-edge_uni[3]*hvy)* gam/(4*pi)
            v_dou[2] += temp * (-edge_uni[1]*hvz+edge_uni[3]*hvx)* gam/(4*pi)
            v_dou[3] += temp * (edge_uni[1]*hvy-edge_uni[2]*hvx)* gam/(4*pi)
        else
            #vdou = [0.0;0.0;0.0]
        end
    end
    return
end

function vel_vring!(v_dou::AbstractVector{T}, loc::AbstractVector{T}, gam::T, pts::AbstractMatrix{T}, pans::AbstractVector{U}, 
    npt::U, edge_len::AbstractVector{T}, edge_uni::AbstractMatrix{T},
    r_rankine::Real, r_cutoff::Real) where {T<:Real, U<:Integer}
    vel_line!(v_dou, pts, pans[1:2], edge_uni[:,1], edge_len[1], loc, gam, r_rankine, r_cutoff)
    vel_line!(v_dou, pts, pans[2:3], edge_uni[:,2], edge_len[2], loc, gam, r_rankine, r_cutoff)
    vel_line!(v_dou, pts, pans[3:4], edge_uni[:,3], edge_len[3], loc, gam, r_rankine, r_cutoff)
    vel_line!(v_dou, pts, [pans[4];pans[1]], edge_uni[:,4], edge_len[4], loc, gam, r_rankine, r_cutoff)
    return
end

function vel_vring2!(v_dou::AbstractVector{T}, loc::AbstractVector{T}, gam::T, pts::AbstractMatrix{T}, pans::AbstractVector{U}, 
    npt::U, edge_len::AbstractVector{T}, edge_uni::AbstractMatrix{T},
    r_rankine::Real, r_cutoff::Real) where {T<:Real, U<:Integer}
    vel_line!(v_dou, pts, pans[1:2], edge_uni[:,1], edge_len[1], loc, gam, r_rankine, r_cutoff)
    vel_line!(v_dou, pts, pans[2:3], edge_uni[:,2], edge_len[2], loc, gam, r_rankine, r_cutoff)
    vel_line!(v_dou, pts, pans[3:4], edge_uni[:,3], edge_len[3], loc, gam, r_rankine, r_cutoff)
    #vel_line!(v_dou, pts, [pans[4];pans[1]], edge_uni[:,4], edge_len[4], loc, gam, r_rankine, r_cutoff)
    return
end

function vel_vring3!(v_dou::AbstractVector{T}, loc::AbstractVector{T}, gam::T, pts::AbstractMatrix{T}, pans::AbstractVector{U}, 
    npt::U, edge_len::AbstractVector{T}, edge_uni::AbstractMatrix{T},
    r_rankine::Real, r_cutoff::Real) where {T<:Real, U<:Integer}
    vel_line!(v_dou, pts, pans[1:2], edge_uni[:,1], edge_len[1], loc, gam, r_rankine, r_cutoff)
    #vel_line!(v_dou, pts, pans[2:3], edge_uni[:,2], edge_len[2], loc, gam, r_rankine, r_cutoff)
    vel_line!(v_dou, pts, pans[3:4], edge_uni[:,3], edge_len[3], loc, gam, r_rankine, r_cutoff)
    vel_line!(v_dou, pts, [pans[4];pans[1]], edge_uni[:,4], edge_len[4], loc, gam, r_rankine, r_cutoff)
    return
end

function vel_vring4!(v_dou::AbstractVector{T}, loc::AbstractVector{T}, gam::T, pts::AbstractMatrix{T}, pans::AbstractVector{U}, 
    npt::U, edge_len::AbstractVector{T}, edge_uni::AbstractMatrix{T},
    r_rankine::Real, r_cutoff::Real) where {T<:Real, U<:Integer}
    vel_line!(v_dou, pts, pans[1:2], edge_uni[:,1], edge_len[1], loc, gam, r_rankine, r_cutoff)
    #vel_line!(v_dou, pts, pans[2:3], edge_uni[:,2], edge_len[2], loc, gam, r_rankine, r_cutoff)
    vel_line!(v_dou, pts, pans[3:4], edge_uni[:,3], edge_len[3], loc, gam, r_rankine, r_cutoff)
    #vel_line!(v_dou, pts, [pans[4];pans[1]], edge_uni[:,4], edge_len[4], loc, gam, r_rankine, r_cutoff)
    return
end

function vel_vring5!(v_dou::AbstractVector{T}, loc::AbstractVector{T}, gam::T, pts::AbstractMatrix{T}, pans::AbstractVector{U}, 
    npt::U, edge_len::AbstractVector{T}, edge_uni::AbstractMatrix{T},
    r_rankine::Real, r_cutoff::Real) where {T<:Real, U<:Integer}
    vel_line!(v_dou, pts, pans[1:2], edge_uni[:,1], edge_len[1], loc, gam, r_rankine, r_cutoff)
    vel_line!(v_dou, pts, pans[2:3], edge_uni[:,2], edge_len[2], loc, gam, r_rankine, r_cutoff)
    #vel_line!(v_dou, pts, pans[3:4], edge_uni[:,3], edge_len[3], loc, gam, r_rankine, r_cutoff)
    vel_line!(v_dou, pts, [pans[4];pans[1]], edge_uni[:,4], edge_len[4], loc, gam, r_rankine, r_cutoff)
    return
end

function compute_vel4!(vels, coord, pansys, wakesys, gam_sol, uinf, r_rankine, r_cutoff)
    npan = pansys.npan
    pan_cpt = pansys.pan_cpt
    pan_norm = pansys.pan_norm
    pan_vert = pansys.pan_vert
    pan_con = pansys.pan_con
    pan_edgelen = pansys.pan_edgelen
    pan_edgeuvec = pansys.pan_edgeuvec
    pan_edgevec = pansys.pan_edgevec
    pan_numpt = pansys.pan_numpt
    pan_neigh_idx = pansys.pan_neigh_idx
    pan_area = pansys.pan_area
    nte = wakesys.nte
    te_pan_idx = wakesys.te_pan_idx
    te_con = wakesys.te_con
    te_vert = wakesys.te_vert
    te_numpt = wakesys.te_numpt
    te_edgelen = wakesys.te_edgelen
    te_edgeuvec = wakesys.te_edgeuvec
    airfoil_strips = pansys.airfoil_strips
    for j = 1:npan
        @views vel_vring!(vels, coord, gam_sol[j], pan_vert, pan_con[:,j], pan_numpt[j], pan_edgelen[:,j], pan_edgeuvec[:,:,j], r_rankine, r_cutoff)
    end
    for j = 1:nte
        @views vel_vring!(vels, coord, gam_sol[te_pan_idx[j]], te_vert, te_con[:,j], te_numpt[j], te_edgelen[:,j], te_edgeuvec[:,:,j], r_rankine, r_cutoff)
    end
    vels[:] .+= uinf
    return
end