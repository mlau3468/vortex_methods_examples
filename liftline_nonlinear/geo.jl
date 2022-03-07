using LinearAlgebra
using Statistics

function buildRectHShoe(span, chord, n)
    # defined clockwise
    dy = span/n

    panVerts = zeros(Float64, 3, 2*(n+1))
    panCon = zeros(Int64,4, n)
    panCpts = zeros(Float64, 3, n)# collocation points
    bndLen = zeros(n)
    chordDir = zeros(3,n)

    le_loc = 0.25

    span_mult = 50 # length of trailing element relative to span

    npt = 1
    for i = 1:n
        p1 = [le_loc*chord + chord + span*span_mult; dy*(i-1); 0]
        p2 = [le_loc*chord; dy*(i-1); 0]
        p3 = [le_loc*chord; dy*(i); 0]
        p4 = [le_loc*chord + chord ; dy*(i); 0]

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
            panVerts[:,npt] = p3
            panVerts[:,npt+1] = p4
            
            panCon[1,i] =  panCon[4,i-1]
            panCon[2,i] =  panCon[3,i-1]
            panCon[3,i] = npt
            panCon[4,i] = npt+1

            npt = npt + 2
            
        end
    end

    return panVerts, panCon, panCpts, bndLen, chordDir
end

function calcHshoeNorm(panVerts, panCon)
    npan = size(panCon,2)
    panNorms = zeros(3,npan)
    for i = 1:npan
        v1 = panVerts[:,panCon[3,i]] - panVerts[:,panCon[1,i]]
        v2 = panVerts[:,panCon[2,i]] - panVerts[:,panCon[4,i]]
        newNorm = cross(v1, v2)
        newNorm = newNorm ./ norm(newNorm)
        panNorms[:,i] .= newNorm
    end
    return panNorms
end

function createRect(span, chord, nspan, nchord, offset, pitch)
    # panels are defined counter clockwise
    pitch = deg2rad(pitch)
    pitch_loc = [0.25*chord; 0; 0]
    rot = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)]
    # create geometry
    te_idx = Int64[]
    pts = zeros(Float64, 3,4*nchord*nspan)
    panels = zeros(Int64,4,nchord*nspan)
    npt = 0 # points index counter
    npan = 0 # panels index counter
    for i = 0:nchord-1
        for j = 0:nspan-1
            p1 = [i*chord/nchord; j*span/nspan; 0]
            p2 = [(i+1)*chord/nchord; j*span/nspan; 0]
            p3 = [(i+1)*chord/nchord; (j+1)*span/nspan; 0]
            p4 = [i*chord/nchord; (j+1)*span/nspan; 0]
            
            
            pts_new = [p1 p2 p3 p4]
            # apply rotation
            for k=1:4
                p = pts_new[:,k]
                p = p .- pitch_loc
                p = rot * p
                p = p + pitch_loc
                pts_new[:,k] = p .+ offset
            end

            # store data
            npan = npan + 1
            for k = 1:4
                npt = npt + 1
                pts[:,npt] = pts_new[:,k]
                panels[k, npan] = npt
            end
            if i+1==nchord
                push!(te_idx, i*nspan + j + 1)
            end
        end
    end
    # merge vertices
    pts, panels = mergeVerts(pts, panels, 1e-6)

    return pts, panels
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


function readCart3D(trifile, scale=1)
    lines = readlines(trifile)
    line1 = split(lines[1], " ", keepempty=false)
    npan = parse(Int64, line1[2])
    nvert = parse(Int64, line1[1])

    pts = zeros(Float64,3,nvert)
    for i = 1:nvert
        pt = split(lines[1+i], " ", keepempty=false)
        pts[:,i] = parse.(Float64, pt)
    end

    pans = zeros(Int64,4,npan)
    for i = 1:npan
        pt = split(lines[1+i+nvert], " ", keepempty=false)
        pans[1:3,i] = parse.(Int64, pt)
    end
    pts = pts .* scale
    println("Number of Panels: $npan")
    println("Number of Vertices: $nvert")
    return pts, pans, npan
end

function readDegenGeo(fname, tol=1e-6, scale=1)
    lines = readlines(fname)
    i = 1
    verts = []
    pans = []
    next_qua = (2,3,4,1)
    # build geometry
    skippt = 0
    while i < length(lines) + 1
        if occursin("SURFACE_NODE" , lines[i])
            flipnorm = cmp(split(lines[i-2], ",")[7], "0")
            temp = split(lines[i], ",")
            nxsec = parse(Int64,temp[2])
            ptsec = parse(Int64,temp[3])
            i = i + 2
            # get vertices
            for j = 1:nxsec*ptsec
                temp = split(lines[i], ",")
                new_pt = [parse(Float64, temp[1]);parse(Float64, temp[2]);parse(Float64, temp[3])]
                if length(verts) == 0
                    verts = new_pt
                else
                    verts = hcat(verts, new_pt)
                end
                i = i + 1
            end
            
            # get panels
            for n = 1:nxsec-1
                for j = 1:ptsec-1
                    new_pan = [(n-1)*ptsec+j; n*ptsec + j; n*ptsec + j+1; (n-1)*ptsec + j+1]
                    if flipnorm == 1
                        new_pan = [new_pan[1];new_pan[4];new_pan[3];new_pan[2]]
                    end
                    new_pan = new_pan .+ skippt

                    # if some points repeated, this is a triangle
                    x = 0
                    new_pan2 = [0;0;0;0]
                    for k = 1:4
                        if ~samePt(verts[:,new_pan[k]], verts[:,new_pan[next_qua[k]]],tol)
                            x = x + 1
                            new_pan2[x] = new_pan[k]
                        end
                    end

                    # append to panel list
                    if length(pans) == 0
                        pans = new_pan2
                    else
                        pans = hcat(pans, new_pan2)
                    end
                end
            end
            skippt = size(verts, 2)
        else
            i = i + 1
        end
    end
    verts, pans = mergeVerts(verts, pans, tol)
    verts = verts.*scale
    verts = convert.(Float32, verts)
    pans = convert.(Int32, pans)

    npan = size(pans, 2)
    nvert = size(verts, 2)
    println("Number of Panels: $npan")
    println("Number of Vertices: $nvert")
    return verts, pans, npan
end

function ptinList(ptlist, newpt, tol=1e-6)
    idx = 1
    found = false
    while !found && idx<=size(ptlist, 2)
        if samePt(ptlist[:,idx], newpt, tol)
            found = true
        else
            idx = idx + 1
        end    
    end
    if !found
       return 0
    else
        return idx
    end
end

function mergeVerts(verts, pans, tol=1e-6)
    new_verts = []
    new_pans = []
    for i = 1:size(pans,2)
        # check if this panel is a triangle
        if pans[4,i] == 0
            nvert = 3
            pts = verts[:,pans[1:3,i]]
        else
            nvert = 4
            pts = verts[:,pans[:,i]]
        end

        if length(new_verts) == 0
            if nvert == 4
                new_verts = verts[:,pans[:,i]]
                new_pans = [1;2;3;4]
            elseif nvert == 3
                new_verts = verts[:,pans[1:3,i]]
                new_pans = [1;2;3;0]
            end
        else
            idx = [0;0;0;0]
            for j = 1:nvert
                res = ptinList(new_verts, pts[:,j], tol)
                if res == 0
                    new_verts = hcat(new_verts, pts[:,j])
                    idx[j] = size(new_verts,2)
                else
                    idx[j] = res
                end
            end
            new_pans = hcat(new_pans, idx)
        end
    end
    nmerge = size(verts,2) - size(new_verts, 2)
    println("Merged $nmerge points")
    return new_verts, new_pans
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

function calcPanProps!(panVert, panCon, panCpt, panNorm, panNPt, 
    panEdgeVec, panEdgeLen, panEdgeUVec, panArea, panTang, panSinTi, panCosTi)

    next_tri = [2 3 1]
    next_qua = [2 3 4 1]

    for i = 1:size(panCon,2)

        # check if panel has 3 or 4 points
        if panCon[4,i] == 0
            panNPt[i] = 3
        else
            panNPt[i] = 4
        end

        if panNPt[i] == 3
            panCpt[:,i] = mean(panVert[:,panCon[1:3,i]],dims=2)
        else
            panCpt[:,i] = mean(panVert[:,panCon[:,i]],dims=2)
        end

        # edge vectors
        if panNPt[i] == 3
            for j = 1:3
                panEdgeVec[:,j,i] = panVert[:, panCon[next_tri[j], i]] .- panVert[:,panCon[j,i]]
            end
        else
            for j = 1:4
                panEdgeVec[:,j,i] = panVert[:, panCon[next_qua[j], i]] .- panVert[:,panCon[j,i]]
            end
        end

        # edge lengths
        for j = 1:panNPt[i]
            panEdgeLen[j,i] = norm(panEdgeVec[:,j,i])
            panEdgeUVec[:,j,i] = panEdgeVec[:,j,i] ./ panEdgeLen[j,i]
        end

        # normals
        v3 = cross(panEdgeVec[:,1,i], panEdgeVec[:,2,i])
        panNorm[:,i] = v3./norm(v3)

        # area
        if panNPt[i] == 3
            panArea[i] = 0.5*norm(v3)
        else
            panArea[i] = norm(v3)
        end

        # local tangent unit vector as in PANAIR
        tanl = 0.5 .* (panVert[:,panCon[panNPt[i],i]] .+ panVert[:,panCon[1,i]]) .- panCpt[:,i]
        panTang[:,1,i] = tanl ./ norm(tanl)
        panTang[:,2,i] = cross(panNorm[:,i], panTang[:,1,i])

        #sinti and costi
        for j = 1:panNPt[i]
            panCosTi[j,i] = sum(panEdgeUVec[:,j,i] .* panTang[:,1,i])
            panSinTi[j,i] = sum(panEdgeUVec[:,j,i] .* panTang[:,2,i])
        end
    end    
    return
end

function calcPanProps(panVert, panCon)
    npan = size(panCon, 2)
    panCpt = zeros(Float32, 3, npan)
    panNorm = zeros(Float32, 3, npan)
    panNPt = zeros(Int32, npan)
    panEdgeVec = zeros(Float32, 3, 4, npan)
    panEdgeLen = zeros(Float32, 4, npan)
    panEdgeUVec = zeros(Float32, 3, 4, npan)
    panTang = zeros(Float32, 3, 2, npan)
    panCosTi = zeros(Float32, 4, npan)
    panSinTi = zeros(Float32, 4, npan)
    panArea = zeros(Float32, npan)
    calcPanProps!(panVert, panCon, panCpt, panNorm, panNPt, panEdgeVec, panEdgeLen, panEdgeUVec, panArea, panTang, panSinTi, panCosTi)
    return panCpt, panNorm, panNPt, panEdgeVec, panEdgeLen, panEdgeUVec, panArea, panTang, panSinTi, panCosTi
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

function calcLineProp!(linePos, lineEdgeVec, lineEdgeLen, lineEdgeUVec)
    for i = 1:size(linePos,3)
        #lineEdgeVec[:,i] = linePos[:,2,i] .- linePos[:,1,i]
        lineEdgeVec[1,i] = linePos[1,2,i] - linePos[1,1,i]
        lineEdgeVec[2,i] = linePos[2,2,i] - linePos[2,1,i]
        lineEdgeVec[3,i] = linePos[3,2,i] - linePos[3,1,i]

        #lineEdgeLen[i] = norm(lineEdgeVec[:,i])
        lineEdgeLen[i] = sqrt(lineEdgeVec[1,i]^2 + lineEdgeVec[2,i]^2 + lineEdgeVec[3,i]^2)
        
        #lineEdgeUVec[:,i] = lineEdgeVec[:,i] ./ lineEdgeLen[i]
        lineEdgeUVec[1,i] = lineEdgeVec[1,i] ./ lineEdgeLen[i]
        lineEdgeUVec[2,i] = lineEdgeVec[2,i] ./ lineEdgeLen[i]
        lineEdgeUVec[3,i] = lineEdgeVec[3,i] ./ lineEdgeLen[i]
    end
    return
end

function extendPtList(ptlist, newpt, tol=1e-8)
    # if newpt is in list, returns the index in the list
    # if newpt not in list, adds to list, then returns the index of the new point
    res = ptinList(ptlist, newpt, tol)
    if res == 0 # not in point list
        ptlist = hcat(ptlist, newpt)
        idx = size(ptlist,2)
    else # point in list
        idx = res
    end
    return ptlist, idx
end

function createTE!(verts, pans, pan_norms, pan_npt, pan_neighIdx, pan_neighSide, pan_neighDir, pan_neighN, pan_cpt, uinf, dt, te_scale, dirvec)
    # find trailing edges
    nte = 0
    te_angle = 30
    npan = size(pans, 2)
    flags = zeros(4, npan)
    prev_qua = [4 1 2 3]
    next_qua = [2 3 4 1]
    prev_tri = [3 1 2]
    next_tri = [2 3 1]

    teVert = []
    teCon = []
    wake_panidx = []
    wake_dir = []
    for i = 1:npan
        for e1 = 1:pan_npt[i]
            i2 = pan_neighIdx[e1,i]
            if i2 > 0  && flags[e1,i] == 0# if has neighbor and no TE created from this edge
                n1 = pan_norms[:,i]
                n2 = pan_norms[:,i2]
                
                temp = dot(n1,n2)
                if temp < 0 && temp > -1
                    a = 180-rad2deg(acos(temp))
                    # check for convexity
                    v1 = pan_cpt[:,i2] .- pan_cpt[:,i] # vector from 1st panel to 2nd panel center point
                    convex = dot(v1, n1) < 0
                    if a <= te_angle && convex # this is a trailing edge
                        e2 = pan_neighSide[e1,i]
                        flags[e1,i] = 1
                        flags[e2,i2] = 1
                        
                        if pan_npt[i] == 4
                            nextPt = next_qua
                            prevPt = prev_qua
                        else
                            nextPt = next_tri
                            prevPt = prev_tri
                        end

                        # remove neighbor
                        pan_neighIdx[e1,i] = 0
                        pan_neighSide[e1,i] = 0
                        pan_neighDir[e1,i] = 0
                        pan_neighN[i] = pan_neighN[i] - 1

                        pan_neighIdx[e2,i2] = 0
                        pan_neighSide[e2,i2] = 0
                        pan_neighDir[e2,i2] = 0
                        pan_neighN[i2] = pan_neighN[i2] - 1

                        p1 = verts[:,pans[nextPt[e1],i]]
                        p2 = verts[:,pans[e1,i]]
                        vec = norm(uinf).*dt.*te_scale.*dirvec
                        
                        p3 = p2 .+ vec
                        p4 = p1 .+ vec
                        
                        pts = [p1 p2 p3 p4]
                        
                        # add wake vertices and panel
                        if length(teVert) == 0
                            teVert = pts
                            idx = [1;2;3;4]
                        else
                            idx = [0;0;0;0]
                            for k = 1:4
                                teVert, idx[k] = extendPtList(teVert, pts[:,k])
                            end
                        end

                        if length(teCon) == 0
                            teCon = idx
                            wake_dir = vec
                            wake_panidx = [i;i2]
                        else
                            teCon = hcat(teCon, idx)
                            wake_dir = hcat(wake_dir, vec)
                            wake_panidx = hcat(wake_panidx, [i;i2])
                        end
                        nte = nte + 1
                    end
                end
            end
        end
    end
    println("Number of Trailing Edges: $nte")

    convert.(Float32, teVert)
    convert.(Int32, teCon)
    convert.(Int32, wake_panidx)
    convert.(Float32, wake_dir)

    tecoords1 = zeros(Float32, 3, nte)
    tecoords2 = zeros(Float32, 3, nte)

    for i = 1:nte
        tecoords1[:,i] .= teVert[:,teCon[3,i]]
        tecoords2[:,i] .= teVert[:,teCon[4,i]]
    end

    return teVert, teCon, wake_panidx, wake_dir, tecoords1, tecoords2, nte
end

function calc_CHTLS_stencil!(pans, cpt, norms,neighN, neighIdx, npts, panCHTLS)
    for i = 1:size(pans,2)
        n_neigh = neighN[i]
        #A differences, B constraints, W weights
        A = zeros(n_neigh, 3)
        W = zeros(n_neigh+1, n_neigh+1)
        B = zeros(1,3)
        B[1,:] = norms[:,i]
        i_n = 0

        for i_nn = 1:npts[i]
            if neighIdx[i_nn,i] > 0
                i_n = i_n + 1
                dx =  cpt[:,neighIdx[i_nn,i]] - cpt[:,i]
                A[i_n,:] = dx
                W[i_n,i_n] = 1/norm(dx)*max(0, abs(sum(norms[:,i].*norms[:,neighIdx[i_nn,i]])))
            end
        end

        W[n_neigh+1, n_neigh+1] = sum(W)./n_neigh
        V = zeros(3,1)
        V[:,1] = B[1,:]
        Q, R = qr(V)
        C = zeros(n_neigh+1, 3)
        C[1:n_neigh,:] = A
        C[n_neigh+1,:] = B[1,:]
        CQ = C*Q
        Cls_tilde = transpose(CQ[:,2:3]) * (W*CQ[:,2:3])
        # take inverse
        iCls_tilde = inv(Cls_tilde)

        r1 = R[1,1]
        chtls_temp = zeros(3, n_neigh+1)
        chtls_temp[1,n_neigh+1] = 1/R[1,1]
        chtls_temp[2:3,1:n_neigh] = iCls_tilde * (transpose(CQ[1:n_neigh,2:3]) * W[1:n_neigh, 1:n_neigh])
        chtls_temp[2:3, n_neigh+1:n_neigh+1] = iCls_tilde * (
        (transpose(CQ[n_neigh+1:n_neigh+1,2:3]) * W[n_neigh+1:n_neigh+1, n_neigh+1:n_neigh+1]) -
        (transpose(CQ[1:n_neigh+1,2:3]) * W[1:n_neigh+1, 1:n_neigh+1]) * CQ[1:n_neigh+1,1:1]./r1)
        # stencil is defined in local frame
        #chtls_stencil =  R_g' * chtls_stencil 
        panCHTLS[:,1:n_neigh+1,i] = Q * chtls_temp
    end
    return 
end

function mirror_point(pt, planeNorm, planePt)
    k = planeNorm[1]*(-pt[1]+planePt[1]) + planeNorm[2]*(-pt[2]+planePt[2]) + planeNorm[3]*(-pt[3]+planePt[3])
    x3 = 2*planeNorm[1]*k + pt[1]
    y3 = 2*planeNorm[2]*k + pt[2]
    z3 = 2*planeNorm[3]*k + pt[3]
    return [x3;y3;z3]
end

function mirror_point!(pt, planeNorm, planePt, new_pt)
    k = planeNorm[1]*(-pt[1]+planePt[1]) + planeNorm[2]*(-pt[2]+planePt[2]) + planeNorm[3]*(-pt[3]+planePt[3])
    new_pt[1] = 2*planeNorm[1]*k + pt[1]
    new_pt[2] = 2*planeNorm[2]*k + pt[2]
    new_pt[3] = 2*planeNorm[3]*k + pt[3]
    return
end

function reversePans!(pans)
    for i = 1:size(pans, 2)
        if pans[4,i] == 0
            pans[1:3,i] = pans[3:-1:1,i]
        else
            pans[:,i] = pans[4:-1:1,i]
        end
    end
end

function mirror_pans(panVert, panCon, planeNorm, planePt)
    panVertMir = zeros(Float32, 3, size(panVert,2))
    for i = 1:size(panVert,2)
        panVertMir[:,i] = mirror_point(panVert[:,i], planeNorm, planePt)
    end
    panConMir = copy(panCon)
    reversePans!(panConMir)
    return panVertMir, panConMir
end

function mirror_pans!(panVert, panCon, planeNorm, planePt, panVertMir, panConMir)
    for i = 1:size(panVert,2)
        @views mirror_point!(panVert[:,i], planeNorm, planePt, panVertMir[:,i])
    end
    panConMir[:,:] .= panCon[:,:]
    reversePans!(panConMir)
    return
end

function updatePartMir!(partPos, partPosMir, partDir, partDirMir, npart, planeNorm, planePt)
    for i = 1:npart
        partPosMir[:,i] = mirror_point(partPos[:,i], planeNorm, planePt)
        partDirMir[:,i] = -partDir[:,i]
    end
end

function updateLineMir!(linePos, linePosMir, lineEdgeVec, lineEdgeVecMir, lineEdgeUVec, lineEdgeUVecMir, nline, planeNorm, planePt)
    for i = 1:nline
        linePosMir[:,1,i] = mirror_point(linePos[:,2,i], planeNorm, planePt)
        linePosMir[:,2,i] = mirror_point(linePos[:,1,i], planeNorm, planePt)
        lineEdgeVecMir[:,i] = -lineEdgeVec[:,i]
        lineEdgeUVecMir[:,i] = -lineEdgeUVec[:,i]
    end
end

function reflectionMatrix(planeNorm, planePt)
    # Matrix to reflect a point across arbitrary plane

end