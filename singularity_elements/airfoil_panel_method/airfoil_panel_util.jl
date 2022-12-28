# utility functions for airfoil panel methods

function calc_panel_2d_collocation_points(pan_vert::Matrix{<:Real})
    npan = size(pan_vert,2) - 1
    pan_cpt = zeros(2, npan)
    for i = 1:npan
        pan_cpt[:,i] = midpoint2d(pan_vert[:,i], pan_vert[:,i+1])
    end
    return pan_cpt
end

function calc_panel_2d_tangents(pan_vert::Matrix{<:Real})
    # assume panels are counterclockwise, tangent points in ccw direction of indexing
    npan = size(pan_vert,2) - 1
    pan_tan = zeros(2, npan)
    for i = 1:npan
        pan_tan[:,i] = calc_line_tan_2d(pan_vert[:,i], pan_vert[:,i+1])
    end
    return pan_tan
end

function calc_panel_2d_normals(pan_vert::Matrix{<:Real})
    # assume panels are counterclockwise, normal points outwards of airfoil boundary
    npan = size(pan_vert,2) - 1
    pan_norm = zeros(2, npan)
    for i = 1:npan
        pan_norm[:,i] = calc_line_norm_2d(pan_vert[:,i+1], pan_vert[:,i])
    end
    return pan_norm
end

function calc_panel_2d_lengths(pan_vert::Matrix{<:Real})
    npan = size(pan_vert,2) - 1
    pan_len = zeros(npan)
    for i = 1:npan
        pan_len[i] = dist2D(pan_vert[:,i], pan_vert[:,i+1])
    end
    return pan_len
end