
@kwdef struct TeSys{T<:Real, U<:Integer}
    # Panel Mesh
    te_vert::Matrix{T}
    te_con::Matrix{U}

    # Panel Properties
    te_cpt::Matrix{T}
    te_numpt::Vector{U}
    te_edgevec::Array{T,3}
    te_edgelen::Matrix{T}
    te_edgeuvec::Array{T,3}
    nte::U
    ntept::U
    ntevert::U
    te_gam::Vector{T}

    # Panel Neighbors
    te_neigh_idx::Matrix{U}
    te_neigh_side::Matrix{U}
    te_neigh_dir::Matrix{U}
    te_neigh_num::Vector{U}

    # Attached Panel Information
    te_pan_idx::Vector{U}
end

function create_te_vortlat!(pansys, uinf, wakedist, verbose::Bool=false)
    pan_con = pansys.pan_con
    pan_vert = pansys.pan_vert

    # get trailing edge. Note that vertex numering starts at the TE
    npan = size(pan_con, 2)
    te_idx_vlat = 2 # in vortex lattice edge numbering, which edge corresponds to te
    next_qua = SVector{4,Int64}(2,3,4,1)
    te_pan_idx = Int64[]
    te_edge_idx = Int64[]
    for i = 1:npan
        if pansys.pan_neigh_idx[te_idx_vlat,i] == 0
            idx2 = pan_con[te_idx_vlat,i]
            idx1 = pan_con[next_qua[te_idx_vlat],i]
            idx = [idx1;idx2]
            if length(te_edge_idx) == 0
                te_edge_idx = idx
            else
                te_edge_idx = hcat(te_edge_idx, idx)
            end
            push!(te_pan_idx, i)
        end
    end

    # unique te points
    te_pts = unique(te_edge_idx) # downstream points

    # te size
    ntept = length(te_pts)
    nte = size(te_edge_idx,2)
    
    # allocate arrays
    T = eltype(pan_vert)
    te_vert = zeros(T, 3, ntept*2)
    te_con = zeros(Int64, 4, nte)

    for i = 1:ntept
        vec = uinf./norm(uinf) .* wakedist
        te_vert[:,i+ntept] = pan_vert[:,te_pts[i]] # upstream row
        te_vert[:,i] = te_vert[:,i+ntept] + vec # downstream row
    end

    for i = 1:nte
        # build panel connectivity
        idx1 = findfirst(te_edge_idx[1,i] .== te_pts)
        idx2 = findfirst(te_edge_idx[2,i] .== te_pts)
        te_con[1,i] = idx1 + ntept
        te_con[2,i] = idx2 + ntept
        te_con[3,i] = idx2 
        te_con[4,i] = idx1
    end

    # calculate wake panel properties
    te_cpt, te_numpt, te_edgevec, te_edgelen, te_edgeuvec = calc_wakepanprop(te_vert, te_con)
    te_neigh_idx, te_neigh_side, te_neigh_dir, te_neigh_num = calc_neighbors(te_con, te_numpt)

    if verbose
        println("Number of Trailing Edges: $nte")
    end

    ntevert = size(te_vert, 2)

    return TeSys(
        # Panel Mesh
        te_vert = te_vert,
        te_con = te_con,
    
        # Panel Properties
        te_cpt = te_cpt,
        te_numpt = te_numpt,
        te_edgevec = te_edgevec,
        te_edgelen = te_edgelen,
        te_edgeuvec = te_edgeuvec,
        nte = nte,
        ntept = ntept,
        ntevert = ntevert,
        te_gam = zeros(T, nte),
    
        # Panel Neighbors
        te_neigh_idx = te_neigh_idx,
        te_neigh_side = te_neigh_side,
        te_neigh_dir = te_neigh_dir,
        te_neigh_num = te_neigh_num,
    
        # Attached Panel Information
        te_pan_idx = te_pan_idx)
end