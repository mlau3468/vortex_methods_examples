
function samePt(pt1, pt2)
    tol = 1e-6
    return norm(pt1.-pt2) < tol
end

function getNeighbors!(panels)
    idx = [1 2 3 4 1]
    for i = 1:length(panels)
        pan1 = panels[i]
        for ii=1:4
            if pan1.neigh[ii] == 0 # if currently no neighbor, check
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
                                panels[i].neigh[ii] = j
                                panels[j].neigh[jj] = i
                                panels[i].neigh_side[ii] = jj
                                panels[j].neigh_side[jj] = ii
                                panels[i].neigh_dir[ii] = -1
                                panels[j].neigh_dir[jj] = -1
                            elseif (samePt(p1,p11) && samePt(p2,p22))
                                # is neighbor
                                panels[i].neigh[ii] = j
                                panels[j].neigh[jj] = i
                                panels[i].neigh_side[ii] = jj
                                panels[j].neigh_side[jj] = ii
                                panels[i].neigh_dir[ii] = 1
                                panels[j].neigh_dir[jj] = 1
                            end
                        end
                    end
                end
            end
        end
    end
end