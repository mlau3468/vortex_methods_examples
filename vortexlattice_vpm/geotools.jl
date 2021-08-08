
function samePt(pt1, pt2)
    tol = 1e-6
    return norm(pt1.-pt2) < tol
end

function calcneighbors(panels)
    idx = [1 2 3 4 1]
    neigh_idx = zeros(Int64,4,length(panels))
    neigh_side = zeros(Int64,4,length(panels))
    neigh_dir = zeros(Int64,4,length(panels))
    for i = 1:length(panels)
        pan1 = panels[i]
        for ii=1:4
            if neigh_idx[ii,i] == 0 # if currently no neighbor, check
                p1 = pan1.pts[:,idx[ii]]
                p2 = pan1.pts[:,idx[ii+1]]
                for j = 1:length(panels)
                    if i != j
                        pan2 = panels[j]
                        for jj = 1:4
                            p11 = pan2.pts[:,idx[jj]]
                            p22 = pan2.pts[:,idx[jj+1]]
                            if (samePt(p1,p22) && samePt(p2,p11))
                                # is neighbor
                                neigh_idx[ii,i] = j
                                neigh_idx[jj,j] = i
                                neigh_side[ii,i] = jj
                                neigh_side[jj,j] = ii
                                neigh_dir[ii,i] = -1
                                neigh_dir[jj,j] = -1
                            elseif (samePt(p1,p11) && samePt(p2,p22))
                                # is neighbor
                                neigh_idx[ii,i] = j
                                neigh_idx[jj,j] = i
                                neigh_side[ii,i] = jj
                                neigh_side[jj,j] = ii
                                neigh_dir[ii,i] = 1
                                neigh_dir[jj,j] = 1
                            end
                        end
                    end
                end
            end
        end
    end
    return neigh_idx, neigh_side, neigh_dir
end