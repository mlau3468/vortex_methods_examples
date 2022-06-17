using LinearAlgebra
using Statistics

function buildRectHShoeCCW(span, chord, n)
    # defined clockwise
    dy = span/n

    panVerts = zeros(Float64, 3, 2*(n+1))
    panCon = zeros(Int64,4, n)
    panCpts = zeros(Float64, 3, n)# collocation points
    bndLen = zeros(n)
    chordDir = zeros(3,n)

    le_loc = 0.25

    wake_len = 50 # length of wake relative to span

    npt = 0
    for i = 1:n
        p4 = [le_loc*chord + chord + span*wake_len; dy*(i-1); 0]
        p3 = [le_loc*chord; dy*(i-1); 0]
        p2 = [le_loc*chord; dy*(i); 0]
        p1 = [le_loc*chord + chord + span*wake_len; dy*(i); 0]

        # bound vortex segment length
        bndLen[i] = norm(p3-p2)

        # direction vector of chord from LE to TE
        cdir = (p1-p2 + p4-p3)./2
        cdir = cdir ./ norm(cdir)
        chordDir[:,i] .= cdir
         
        # collocation point
        # panel te for collcation point calculation
        p4c = [(le_loc+1)*chord; dy*(i-1); 0]
        p1c = [(le_loc+1); dy*(i); 0]
        panCpts[:,i] = (p1c.+p2.+p3.+p4c)/4
        if i == 1
            panVerts[:,1] = p1
            panVerts[:,2] = p2
            panVerts[:,3] = p3
            panVerts[:,4] = p4
            panCon[:,i] .= [1;2;3;4]
            npt = npt + 4
        else
            panVerts[:,npt+1] = p1
            panVerts[:,npt+2] = p2
            
            panCon[4,i] = panCon[1,i-1]
            panCon[3,i] = panCon[2,i-1]
            panCon[2,i] = npt + 2
            panCon[1,i] = npt + 1

            npt = npt + 2
            
        end
    end
    panNorms = zeros(3,n)
    for i = 1:n
        v1 = panVerts[:,panCon[2,i]] - panVerts[:,panCon[4,i]]
        v2 = panVerts[:,panCon[3,i]] - panVerts[:,panCon[1,i]]
        newNorm = cross(v1, v2)
        newNorm = newNorm ./ norm(newNorm)
        panNorms[:,i] .= newNorm
    end

    return panVerts, panCon, panCpts, bndLen, chordDir, panNorms
end

function buildRectHShoe(span, chord, n)
    # defined clockwise
    dy = span/n

    panVerts = zeros(Float64, 3, 2*(n+1))
    panCon = zeros(Int64,4, n)
    panCpts = zeros(Float64, 3, n)# collocation points
    bndLen = zeros(n)
    chordDir = zeros(3,n)

    le_loc = 0.25

    wake_len = 50 # length of wake relative to span

    npt = 0
    for i = 1:n
        p1 = [le_loc*chord + chord + span*wake_len; dy*(i-1); 0]
        p2 = [le_loc*chord; dy*(i-1); 0]
        p3 = [le_loc*chord; dy*(i); 0]
        p4 = [le_loc*chord + chord + span*wake_len; dy*(i); 0]

        # panel te for collcation point calculation
        p1c = [(le_loc+1)*chord; dy*(i-1); 0]
        p4c = [(le_loc+1); dy*(i); 0]

        # bound vortex segment length
        bndLen[i] = norm(p3-p2)

        # direction vector of chord from LE to TE
        cdir = (p1-p2 + p4-p3)./2
        cdir = cdir ./ norm(cdir)
        chordDir[:,i] .= cdir
         

        # collocation point
        panCpts[:,i] = (p1c.+p2.+p3.+p4c)/4
        if i == 1
            panVerts[:,1] = p1
            panVerts[:,2] = p2
            panVerts[:,3] = p3
            panVerts[:,4] = p4
            panCon[:,i] .= [1;2;3;4]
            npt = npt + 4
        else
            panVerts[:,npt+1] = p3
            panVerts[:,npt+2] = p4
            
            panCon[1,i] =  panCon[4,i-1]
            panCon[2,i] =  panCon[3,i-1]
            panCon[3,i] = npt + 1
            panCon[4,i] = npt + 2

            npt = npt + 2
            
        end
    end
    panNorms = zeros(3,n)
    for i = 1:n
        v1 = panVerts[:,panCon[3,i]] - panVerts[:,panCon[1,i]]
        v2 = panVerts[:,panCon[2,i]] - panVerts[:,panCon[4,i]]
        newNorm = cross(v1, v2)
        newNorm = newNorm ./ norm(newNorm)
        panNorms[:,i] .= newNorm
    end

    return panVerts, panCon, panCpts, bndLen, chordDir, panNorms
end

function buildRect(span, chord, n)
    # defined clockwise
    dy = span/n

    panVerts = zeros(Float64, 3, 2*(n+1))
    panCon = zeros(Int64,4, n)
    panCpts = zeros(Float64, 3, n)# collocation points
    bndLen = zeros(n)
    chordDir = zeros(3,n)
    panArea = zeros(n)
    panNorms = zeros(3,n)

    le_loc = 0.25

    extend_len = 0 # extend vortex line downstream by this amount, =0 for vortex ring

    npt = 0
    for i = 1:n
        # use this to extend panel way downstream to mimick horshoe
        p4 = [le_loc*chord + chord + span*extend_len; dy*(i-1); 0]
        p3 = [le_loc*chord; dy*(i-1); 0]
        p2 = [le_loc*chord; dy*(i); 0]
        p1 = [le_loc*chord + chord + span*extend_len; dy*(i); 0]

        # use this for actual vortex ring
        p4 = [le_loc*chord + chord; dy*(i-1); 0]
        p3 = [le_loc*chord; dy*(i-1); 0]
        p2 = [le_loc*chord; dy*(i); 0]
        p1 = [le_loc*chord + chord; dy*(i); 0]

        # bound vortex segment length
        bndLen[i] = norm(p3-p2)

        # direction vector of chord from LE to TE
        cdir = (p1-p2 + p4-p3)./2
        cdir = cdir ./ norm(cdir)
        chordDir[:,i] .= cdir
         
        # collocation point
        # panel te for collcation point calculation
        p4c = [(le_loc+1)*chord; dy*(i-1); 0]
        p1c = [(le_loc+1); dy*(i); 0]
        panCpts[:,i] = (p1c.+p2.+p3.+p4c)/4
        if i == 1
            panVerts[:,1] = p1
            panVerts[:,2] = p2
            panVerts[:,3] = p3
            panVerts[:,4] = p4
            panCon[:,i] .= [1;2;3;4]
            npt = npt + 4
        else
            panVerts[:,npt+1] = p1
            panVerts[:,npt+2] = p2
            
            panCon[4,i] = panCon[1,i-1]
            panCon[3,i] = panCon[2,i-1]
            panCon[2,i] = npt + 2
            panCon[1,i] = npt + 1

            npt = npt + 2
            
        end

        # normals
        v1 = p2 .- p4c
        v2 = p3 .- p1c
        v3 = cross(v1, v2)
        panNorms[:,i] .= v3 ./ norm(v3)

        # areas
        panArea[i] = norm(v3)/2
    end

    return panVerts, panCon, panCpts, bndLen, chordDir, panNorms, panArea
end

function deg180(a)
    # takes angle a in degrees, converts to be between -180 and 180 degrees
    a = mod(a, 360)
    if a > 180
        a = a - 360
    elseif a < -180
        a = a + 360
    end
    return a
end

# ALL BELOW IS COPIED FROM VPM CODE -----------------------------------


function eulerMat(a, b, g, degrees=false)
    if degrees
        a = deg2rad(a)
        b = deg2rad(b)
        g = deg2rad(g)
    end
    # computes euler matrix to perform 321 rotation
    return [[cos(a)*cos(b) cos(a)*sin(b)*sin(g)-sin(a)*cos(g) cos(a)*sin(b)*cos(g)+sin(a)*sin(g)]; 
    [sin(a)*cos(b) sin(a)*sin(b)*sin(g)+cos(a)*cos(g) sin(a)*sin(b)*cos(g)-cos(a)*sin(g)];
    [-sin(b) cos(b)*sin(g) cos(b)*cos(g)]]  
end


function samePt(pt1, pt2, tol=1e-6)
    return norm(pt1.-pt2) < tol
end

function calcneighbors!(panels, panNPt, neigh_idx, neigh_side, neigh_dir, n_neigh)
    npan = size(panels,2)
    next_qua = [2 3 4 1]
    next_tri = [2 3 1]
    @views for i = 1:npan
        @views for ii=1:panNPt[i]
            if neigh_idx[ii,i] == 0 # if currently no neighbor, check
                p1 = panels[ii,i]
                if panNPt[i] == 4
                    p2 = panels[next_qua[ii],i]
                else
                    p2 = panels[next_tri[ii],i]
                end
                @views for j = 1:npan
                    if i != j
                        if panNPt[j] == 4
                            idx2 = next_qua
                        else
                            idx2 = next_tri
                        end
                        @views for jj = 1:panNPt[j]
                            p11 = panels[jj,j]
                            p22 = panels[idx2[jj],j]
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

function calcneighbors(panCon, panNPt)
    npan = size(panCon, 2)
    panNeighIdx = zeros(Int32, 4, npan)
    panNeighSide = zeros(Int32, 4, npan)
    panNeighDir = zeros(Int32, 4, npan)
    panNNeigh = zeros(Int32, npan)
    calcneighbors!(panCon, panNPt, panNeighIdx, panNeighSide, panNeighDir, panNNeigh)
    return  panNeighIdx, panNeighSide, panNeighDir, panNNeigh
end
