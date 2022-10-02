function calc_neighbors(pan_con::AbstractMatrix{Int64}, pan_numpt::AbstractVector{Int64})
    npan = size(pan_con, 2)
    pan_neigh_idx = zeros(Int64, 4, npan)
    pan_neigh_side = zeros(Int64, 4, npan)
    pan_neigh_dir = zeros(Int64, 4, npan)
    pan_neigh_num = zeros(Int64, npan)
    calc_neighbors!(pan_neigh_idx, pan_neigh_side, pan_neigh_dir, pan_neigh_num, pan_con, pan_numpt, npan)
    return  pan_neigh_idx, pan_neigh_side, pan_neigh_dir, pan_neigh_num
end

function calc_neighbors!(neigh_idx::AbstractMatrix{T}, neigh_side::AbstractMatrix{T}, neigh_dir::AbstractMatrix{T}, n_neigh::AbstractVector{T}, pan_con::AbstractMatrix{T}, pan_numpt::AbstractVector{T}, npan::Integer) where {T<:Integer}
    next_qua = SVector{4,Int64}(2,3,4,1)
    next_tri = SVector{3,Int64}(2,3,1)
    @views @inbounds for i = 1:npan
        @views for ii=1:pan_numpt[i]
            if neigh_idx[ii,i] == 0 # if currently no neighbor, check
                p1 = pan_con[ii,i]
                if pan_numpt[i] == 4
                    p2 = pan_con[next_qua[ii],i]
                else
                    p2 = pan_con[next_tri[ii],i]
                end
                @views for j = 1:npan
                    if i != j
                        if pan_numpt[j] == 4
                            idx2 = next_qua
                        else
                            idx2 = next_tri
                        end
                        @views for jj = 1:pan_numpt[j]
                            p11 = pan_con[jj,j]
                            p22 = pan_con[idx2[jj],j]
                            if p1 == p22 && p2 == p11
                                neigh_idx[ii,i] = j
                                neigh_idx[jj,j] = i
                                neigh_side[ii,i] = jj
                                neigh_side[jj,j] = ii
                                neigh_dir[ii,i] = -1
                                neigh_dir[jj,j] = -1
                                n_neigh[i] = n_neigh[i] + 1
                                n_neigh[j] = n_neigh[j] + 1
                            elseif p1 == p11 && p2 == p22
                                neigh_idx[ii,i] = j
                                neigh_idx[jj,j] = i
                                neigh_side[ii,i] = jj
                                neigh_side[jj,j] = ii
                                neigh_dir[ii,i] = 1
                                neigh_dir[jj,j] = 1
                                n_neigh[i] = n_neigh[i] + 1
                                n_neigh[j] = n_neigh[j] + 1
                            end
                        end
                    end
                end
            end
        end
    end
    return
end

"""
    calc_wakepanprop(pan_vert, pan_con)
Calculates properties for wake panels defined by a collection of vertices panvert and node loops in pan_con.
"""
function calc_wakepanprop(pan_vert::AbstractMatrix{T}, pan_con::AbstractMatrix{Int64}) where {T<:Real}
    npan = size(pan_con, 2)
    pan_numpt = zeros(Int64, npan)
    pan_edgevec = zeros(T, 3, 4, npan)
    pan_edgelen = zeros(T, 4, npan)
    pan_edgeuvec = zeros(T, 3, 4, npan)
    pan_cpt = zeros(T, 3, npan)
    calc_wakepanprop!(pan_cpt, pan_numpt, pan_edgevec, pan_edgelen, pan_edgeuvec, pan_vert, pan_con, npan)
    return pan_cpt, pan_numpt, pan_edgevec, pan_edgelen, pan_edgeuvec
end

function calc_wakepanprop!(wakesys, nwakepan)
    calc_wakepanprop!(wakesys.wake_cpt, wakesys.wake_numpt, wakesys.wake_edgevec, wakesys.wake_edgelen, wakesys.wake_edgeuvec, wakesys.wake_vert, wakesys.wake_con, nwakepan)
    return
end

"""
    calc_wakepanprop!(pan_numpt, pan_edgevec, pan_edgelen, pan_edgeuvec, pan_vert, pan_con, npan)
Calculates properties for wake panels defined by a collection of vertices panvert and node loops in pan_con.
"""
function calc_wakepanprop!(pan_cpt::AbstractMatrix{T}, pan_numpt::AbstractVector{U}, pan_edgevec::AbstractArray{T,3},
                        pan_edgelen::AbstractMatrix{T}, pan_edgeuvec::AbstractArray{T,3}, pan_vert::AbstractMatrix{T}, 
                        pan_con::AbstractMatrix{U}, npan::Integer) where {T<:Real,U<:Integer}
    # panels defined couterclockwise when looking from above
    next_tri = [2;3;1]
    next_qua = [2;3;4;1]

    calc_pannpt!(pan_con, pan_numpt, npan)
    calc_pancpt!(pan_cpt, pan_vert, pan_con, npan)

    @views @inbounds Threads.@threads for i = 1:npan
        # edge vectors
        if pan_numpt[i] == 3
            @inbounds for j = 1:3
                # pan_edgevec[:,j,i] = pan_vert[:, pan_con[next_tri[j], i]] .- pan_vert[:,pan_con[j,i]]
                pan_edgevec[1,j,i] = pan_vert[1, pan_con[next_tri[j], i]] - pan_vert[1,pan_con[j,i]]
                pan_edgevec[2,j,i] = pan_vert[2, pan_con[next_tri[j], i]] - pan_vert[2,pan_con[j,i]]
                pan_edgevec[3,j,i] = pan_vert[3, pan_con[next_tri[j], i]] - pan_vert[3,pan_con[j,i]]

                pan_edgelen[j,i] = norm(pan_edgevec[:,j,i])

                # pan_edgeuvec[:,j,i] = pan_edgevec[:,j,i] ./ pan_edgelen[j,i]
                pan_edgeuvec[1,j,i] = pan_edgevec[1,j,i] / pan_edgelen[j,i]
                pan_edgeuvec[2,j,i] = pan_edgevec[2,j,i] / pan_edgelen[j,i]
                pan_edgeuvec[3,j,i] = pan_edgevec[3,j,i] / pan_edgelen[j,i]
            end
        else
            @inbounds for j = 1:4
                # pan_edgevec[:,j,i] = pan_vert[:, pan_con[next_qua[j], i]] .- pan_vert[:,pan_con[j,i]]
                pan_edgevec[1,j,i] = pan_vert[1, pan_con[next_qua[j], i]] - pan_vert[1,pan_con[j,i]]
                pan_edgevec[2,j,i] = pan_vert[2, pan_con[next_qua[j], i]] - pan_vert[2,pan_con[j,i]]
                pan_edgevec[3,j,i] = pan_vert[3, pan_con[next_qua[j], i]] - pan_vert[3,pan_con[j,i]]

                pan_edgelen[j,i] = norm(pan_edgevec[:,j,i])

                # pan_edgeuvec[:,j,i] = pan_edgevec[:,j,i] ./ pan_edgelen[j,i]
                pan_edgeuvec[1,j,i] = pan_edgevec[1,j,i] / pan_edgelen[j,i]
                pan_edgeuvec[2,j,i] = pan_edgevec[2,j,i] / pan_edgelen[j,i]
                pan_edgeuvec[3,j,i] = pan_edgevec[3,j,i] / pan_edgelen[j,i]
            end
        end
    end    
    return
end

"""
    calc_panprop!(pan_cpt, pan_norm, pan_numpt, pan_edgevec, pan_edgelen, pan_edgeuvec, pan_area, pan_tang, pan_sin_ti, pan_cos_ti, pan_vert, pan_con, npan)
Calculates properties for panels defined by a collection of vertices panvert and node loops in pan_con. Non allocating version.
"""
function calc_panprop!(pan_cpt::AbstractMatrix{T}, pan_norm::AbstractMatrix{T}, pan_numpt::AbstractVector{U}, pan_edgevec::AbstractArray{T,3}, 
                    pan_edgelen::AbstractMatrix{T}, pan_edgeuvec::AbstractArray{T,3}, pan_area::AbstractVector{T}, pan_tang::AbstractArray{T,3}, 
                    pan_sin_ti::AbstractMatrix{T}, pan_cos_ti::AbstractMatrix{T}, pan_vert::AbstractMatrix{T}, pan_con::AbstractMatrix{U}, npan::U) where {T<:Real, U<:Integer}
    # panels defined couterclockwise when looking from above

    next_tri = SVector{3,Int64}(2,3,1)
    next_qua = SVector{4,Int64}(2,3,4,1)

    calc_pancpt!(pan_cpt, pan_vert, pan_con, npan)
    calc_pannpt!(pan_con, pan_numpt, npan)
    @views @inbounds Threads.@threads for i = 1:npan
       
        # edge vectors
        if pan_numpt[i] == 3
            @inbounds for j = 1:3
                pan_edgevec[:,j,i] .= pan_vert[:, pan_con[next_tri[j], i]] .- pan_vert[:,pan_con[j,i]]
            end
        else
            @inbounds for j = 1:4
                pan_edgevec[:,j,i] .= pan_vert[:, pan_con[next_qua[j], i]] .- pan_vert[:,pan_con[j,i]]
            end
        end

        # edge lengths
        @inbounds for j = 1:pan_numpt[i]
            pan_edgelen[j,i] = norm(pan_edgevec[:,j,i])
            pan_edgeuvec[:,j,i] .= pan_edgevec[:,j,i] ./ pan_edgelen[j,i]
        end
        if pan_numpt[i] == 3
            v3 = cross(pan_edgevec[:,1,i], pan_edgevec[:,2,i])
            pan_norm[:,i] .= v3./norm(v3)
            pan_area[i] = 0.5*norm(v3)
        else
            # v1 = pan_vert[:,pan_con[3,i]] .- pan_vert[:,pan_con[1,i]]
            # v2 = pan_vert[:,pan_con[4,i]] .- pan_vert[:,pan_con[2,i]]
            # v3 = cross(v1, v2)
            v1x = pan_vert[1,pan_con[3,i]] - pan_vert[1,pan_con[1,i]]
            v1y = pan_vert[2,pan_con[3,i]] - pan_vert[2,pan_con[1,i]]
            v1z = pan_vert[3,pan_con[3,i]] - pan_vert[3,pan_con[1,i]]
            v2x = pan_vert[1,pan_con[4,i]] - pan_vert[1,pan_con[2,i]]
            v2y = pan_vert[2,pan_con[4,i]] - pan_vert[2,pan_con[2,i]]
            v2z = pan_vert[3,pan_con[4,i]] - pan_vert[3,pan_con[2,i]]
            v3x = v1y*v2z-v1z*v2y
            v3y = v1z*v2x-v1x*v2z
            v3z = v1x*v2y-v1y*v2x
            n = sqrt(v3x^2 + v3y^2 + v3z^2)
            pan_norm[1,i] = v3x/n
            pan_norm[2,i] = v3y/n
            pan_norm[3,i] = v3z/n
            pan_area[i] = 0.5*n
        end
        
        # local tangent unit vector as in PANAIR
        # tanl = 0.5 .* (pan_vert[:,pan_con[pan_numpt[i],i]] .+ pan_vert[:,pan_con[1,i]]) .- pan_cpt[:,i]
        # pan_tang[:,1,i] .= tanl ./ norm(tanl)
        tanlx =  0.5 .* (pan_vert[1,pan_con[pan_numpt[i],i]] .+ pan_vert[1,pan_con[1,i]]) .- pan_cpt[1,i]
        tanly = 0.5 .* (pan_vert[2,pan_con[pan_numpt[i],i]] .+ pan_vert[2,pan_con[1,i]]) .- pan_cpt[2,i]
        tanlz = 0.5 .* (pan_vert[3,pan_con[pan_numpt[i],i]] .+ pan_vert[3,pan_con[1,i]]) .- pan_cpt[3,i]
        tanmag = sqrt(tanlx^2+tanly^2+tanlz^2)
        pan_tang[1,1,i] = tanlx/tanmag
        pan_tang[2,1,i] = tanly/tanmag
        pan_tang[3,1,i] = tanlz/tanmag
        
        # pan_tang[:,2,i] .= cross(pan_norm[:,i], pan_tang[:,1,i])
        pan_tang[1,2,i] = pan_norm[2,i]*pan_tang[3,1,i] - pan_norm[3,i]*pan_tang[2,1,i]
        pan_tang[2,2,i] = pan_norm[3,i]*pan_tang[1,1,i] - pan_norm[1,i]*pan_tang[3,1,i]
        pan_tang[3,2,i] = pan_norm[1,i]*pan_tang[2,1,i] - pan_norm[2,i]*pan_tang[1,1,i]
        
        #sinti and costi
        @inbounds for j = 1:pan_numpt[i]
            #pan_cos_ti[j,i] = sum(pan_edgeuvec[:,j,i] .* pan_tang[:,1,i])
            pan_cos_ti[j,i] = pan_edgeuvec[1,j,i]*pan_tang[1,1,i] + pan_edgeuvec[2,j,i]*pan_tang[2,1,i] + pan_edgeuvec[3,j,i]*pan_tang[3,1,i]
            #pan_sin_ti[j,i] = sum(pan_edgeuvec[:,j,i] .* pan_tang[:,2,i])
            pan_sin_ti[j,i] = pan_edgeuvec[1,j,i]*pan_tang[1,2,i] + pan_edgeuvec[2,j,i]*pan_tang[2,2,i] + pan_edgeuvec[3,j,i]*pan_tang[3,2,i]
        end
    end
    return
end

"""
    calc_panprop(pan_vert, pan_con)
Calculates properties for panels defined by a collection of vertices panvert and node loops in pan_con.
"""
function calc_panprop(pan_vert::AbstractMatrix{T}, pan_con::AbstractMatrix{Int64}) where {T<:Real}
    npan = size(pan_con, 2)
    pan_cpt = zeros(T, 3, npan)
    pan_norm = zeros(T, 3, npan)
    pan_numpt = zeros(Int64, npan)
    pan_edgevec = zeros(T, 3, 4, npan)
    pan_edgelen = zeros(T, 4, npan)
    pan_edgeuvec = zeros(T, 3, 4, npan)
    pan_tang = zeros(T, 3, 2, npan)
    pan_cos_ti = zeros(T, 4, npan)
    pan_sin_ti = zeros(T, 4, npan)
    pan_area = zeros(T, npan)
    calc_panprop!(pan_cpt, pan_norm, pan_numpt, pan_edgevec, pan_edgelen, pan_edgeuvec, pan_area, pan_tang, pan_sin_ti, pan_cos_ti, pan_vert, pan_con, npan)
    return pan_cpt, pan_norm, pan_numpt, pan_edgevec, pan_edgelen, pan_edgeuvec, pan_area, pan_tang, pan_sin_ti, pan_cos_ti
end

"""
    calc_pancpt!(pan_cpt, pan_vert, pan_con, npan)
Calculates panel collocation point.
"""
function calc_pancpt!(pan_cpt::AbstractMatrix{T}, pan_vert::AbstractMatrix{T}, pan_con::AbstractMatrix{Int64}, npan::Int64) where {T<:Real}
    @views @inbounds for i = 1:npan
        if pan_con[4,i] == 0
            #pan_cpt[:,i] = mean(pan_vert[:,pan_con[1:3,i]],dims=2)
            pan_cpt[1,i] = (pan_vert[1,pan_con[1,i]] + pan_vert[1,pan_con[2,i]] + pan_vert[1,pan_con[3,i]])/3
            pan_cpt[2,i] = (pan_vert[2,pan_con[1,i]] + pan_vert[2,pan_con[2,i]] + pan_vert[2,pan_con[3,i]])/3
            pan_cpt[3,i] = (pan_vert[3,pan_con[1,i]] + pan_vert[3,pan_con[2,i]] + pan_vert[3,pan_con[3,i]])/3
        else
            #pan_cpt[:,i] = mean(pan_vert[:,pan_con[:,i]],dims=2)
            pan_cpt[1,i] = (pan_vert[1,pan_con[1,i]] + pan_vert[1,pan_con[2,i]] + pan_vert[1,pan_con[3,i]] + pan_vert[1,pan_con[4,i]])/4
            pan_cpt[2,i] = (pan_vert[2,pan_con[1,i]] + pan_vert[2,pan_con[2,i]] + pan_vert[2,pan_con[3,i]] + pan_vert[2,pan_con[4,i]])/4
            pan_cpt[3,i] = (pan_vert[3,pan_con[1,i]] + pan_vert[3,pan_con[2,i]] + pan_vert[3,pan_con[3,i]] + pan_vert[3,pan_con[4,i]])/4
        end
    end
    return
end

"""
    calc_pannorm!(pan_norm, pan_vert, pan_con, pan_numpt, npan)
Calculates panel unit normal vectors.
"""
function calc_pannorm!(pan_norm::AbstractMatrix{T}, pan_vert::AbstractMatrix{T}, pan_con::AbstractMatrix{Int64}, pan_numpt::AbstractVector{Int64}, npan::Int64) where {T<:Real}
    @views @inbounds for i = 1:npan
        if pan_numpt[i] == 3 # triangle
            v1 = pan_vert[:, pan_con[2,i]] .- pan_vert[:, pan_con[1,i]]
            v2 = pan_vert[:, pan_con[3,i]] .- pan_vert[:, pan_con[2,i]]
            v3 = cross(v1, v2)
            pan_norm[:,i] .= v3./norm(v3)
        else # quad
            v1 = pan_vert[:,pan_con[3,i]] .- pan_vert[:,pan_con[1,i]]
            v2 = pan_vert[:,pan_con[4,i]] .- pan_vert[:,pan_con[2,i]]
            v3 = cross(v1,v2)
            pan_norm[:,i] .= v3./norm(v3)
        end
    end
    return
end

"""
    calc_pannpt!(pan_con, pan_numpt, npan)
Find number of points in each panel.
"""
function calc_pannpt!(pan_con::AbstractMatrix{T}, pan_numpt::AbstractVector{T}, npan::Integer) where {T<:Integer}
    @inbounds for i = 1:npan
        # check if panel has 3 or 4 points
        if pan_con[4,i] == 0
            pan_numpt[i] = 3
        else
            pan_numpt[i] = 4
        end
    end
    return
end