using HDF5
using Statistics
using LinearAlgebra
using DelimitedFiles
using WriteVTK

function factorize2(A)
    # modification of LU decomposition in julia to output single matrix
    L,U = factorize(A)
    res = zeros(size(A))
    for i = 1:size(A,1)
        for j = 1:size(A,1)
            if i <= j
                res[i,j] = U[i,j]
            else
                res[i,j] = L[i,j]
            end
        end
    end
    return res
end

function sourc(pan, loc, dou, rr)
    pts = rr[:, pan.ee]
    
    cpt = pan.center
    radius = norm(loc.-cpt)
    e3 = pan.norm
    zQ = sum((loc.-cpt).*e3)
    Qp = loc .- zQ .* e3

    # get edge lens
    edge_len = pan.edgeLen

    # get edge vectors
    edge_vec = pan.edgeVec

    # settings
    prev_qua = [4 1 2 3]
    next_qua = [2 3 4 1]
    prev_tri = [3 1 2]
    next_tri = [2 3 1]

    farField = false
    if farField
    # far field calculate
    # sou = 
    else
    sou = 0.0
    n_ver = pan.nVert
    for i=1:n_ver
        if n_ver == 3 # triangle element
            indm1 = prev_tri[i]
            indp1 = next_tri[i]
        elseif n_ver == 4 # quad element
            indm1 = prev_qua[i]
            indp1 = next_qua[i]
        end

        R1 = norm(loc .- pts[:,i])
        R2 = norm(loc .- pts[:,indp1])
        den = R1+R2-edge_len[i]


        if den < 1e-6
            println("Too small denominator in source computation with point projection")
            #R1 = norm(loc .- pts[:,1])
            #R2 = norm(loc .- pts[:,indp1])
            #den = R1+R2-edge_len[i]
            
        end
        souLog = log( (R1+R2+edge_len[i])  /   den)
        vi = -sum(cross(Qp.-pts[:,i], edge_vec[:,i]).*e3) / edge_len[i]
        sou = sou + vi*souLog
    end

    sou = sou - zQ * dou
    end
    return sou
end


function dub(pan, loc, rr)
    pts = rr[:, pan.ee]
    cpt = mean(pts, dims=2)
    radius = norm(loc.-cpt)
    e3 = pan.norm

    # get edge vectors
    edge_vec = pan.edgeVec

    farField = false

    if farField
        # far field approximation
    else
        zQ = sum((loc.-cpt).*e3)
        dou = 0.0
        
        n_ver = 4
        for i1 = 1:n_ver
            indp = 1+ (i1 % n_ver)
            indm1 = n_ver - ((n_ver-i1+1) % n_ver)
            ei = - (loc .- pts[:,i1])
            ei = ei ./ norm(ei)
            ap1 = edge_vec[:,i1] .- ei .* sum(ei .* edge_vec[:,i1])
            am1 = -edge_vec[:,indm1] .+ ei .* sum(ei .* edge_vec[:,indm1])
            ap1n = norm(ap1)
            am1n = norm(am1)
            sinB = sum(ei.*cross(am1, ap1)) ./ (ap1n.*am1n)
            cosB = sum(am1.*ap1 ./ (ap1n.*am1n))
            beta = atan(sinB, cosB)
            dou = dou + beta
        end
        #Correct the result to obtain the solid angle (from Gauss-Bonnet theorem)
        if dou < -(n_ver-2) * pi + 1e-5
            dou = dou + (n_ver-2)*pi
        elseif dou > (n_ver-2) * pi - 1e-5
            dou = dou - (n_ver-2)*pi
        end
    end
    return dou
end

function calc_node_vel(r, G,f)
    # calculate velocity of a point whose coordinate is rr. Boundary condition
    #r: point coordinate in global coordinate frame 
    #G: frame rotation rate with respect othe base reference
    #f: frame framve velocity with respect to the base reference
    v = f.+(G*r) # velocity
    return v
end

function vel_dub(pan, loc, rr)
    #=
    !> compute velocity AIC of panel <this>, on the control point <pos>
!! Biot-Savart law. Regular kernel: linear core (Rankine vortex)
!! Compute the velocity induced by a vortex ring (equivalent to a constant
!!  intentsity surface doublet) with intensity 4*pi. <------
=#
    global v_dou = 0.0
    n_ver = pan.nVert
    pts = rr[:, pan.ee]
    # get edge vectors
    edge_vec = pan.edgeVec
    # get edge lens
    edge_len = pan.edgeLen
    # unit vectors
    edge_uni = pan.edgeUni

    r_rankine = 0.1
    r_cutoff = 0.001
    for i = 1:n_ver
        indp1 = 1 + i % n_ver
        indm1 = n_ver - (n_ver-i+1) % n_ver
        av = loc.-pts[:,i]
        ai = sum(av.*edge_uni[:,i])
        r1 = norm(av)
        r2 = norm(loc.-pts[:,indp1])
        hv = av .- ai.*edge_uni[:,i]
        hi = norm(hv)
        
        if hi > edge_len[i]*r_rankine
            global v_dou = v_dou .+ ((edge_len[i].-ai)./r2 .+ ai./r1) ./ (hi.^2) .* cross(edge_uni[:,i], hv)
        else
            if (r1 > edge_len[i]*r_cutoff) && (r2 > edge_len[i]*r_cutoff)
                r_Ran = r_rankine * edge_len[i]
                global v_dou = v_dou .+ ((edge_len[i].-ai)./r2 + ai./r1) ./ (r_Ran.^2) .* cross(edge_uni[:,i], hv)
            end
        end
    end

    return v_dou
end

function vel_sourc(pan, loc, rr)
    phix = 0.0
    phiy = 0.0

    n_ver = pan.nVert
    pts = rr[:,pan.ee]
    # settings
    prev_qua = [4 1 2 3]
    next_qua = [2 3 4 1]
    prev_tri = [3 1 2]
    next_tri = [2 3 1]

    # get edge vectors
    edge_vec = pan.edgeVec
    # get edge lens
    edge_len = pan.edgeLen
    # unit vectors
    edge_uni = pan.edgeUni

    # central point
    cpt = pan.center

    # local tangent unit vector as in PANAIR
    n_sides = pan.nSide
    tang = pan.tang

    sinTi = pan.sinTi
    cosTi = pan.cosTi

    for i = 1:n_ver
        if n_ver == 3 # triangle element
            indm1 = prev_tri[i]
            indp1 = next_tri[i]
        elseif n_ver == 4 # quad element
            indm1 = prev_qua[i]
            indp1 = next_qua[i]
        end

        R1 = norm(loc .- pts[:,i])
        R2 = norm(loc .- pts[:, indp1])

        if (R1+R2-edge_len[i]) < 1e-12
            souLog = 0
        else
            souLog = log((R1+R2+edge_len[i]) / (R1+R2-edge_len[i]))
        end
        phix = phix + sinTi[i] * souLog
        phiy = phiy - cosTi[i] * souLog
    end
    normal = pan.norm
    pdou = dub(pan, loc, rr)

    # missing negative, check equations
    phix = - phix
    phiy = - phiy
    pdou = - pdou

    #println(pdou)
    vel = zeros(3)
    vel[1] = tang[1,1] * phix + tang[1,2]*phiy + normal[1] * pdou
    vel[2] = tang[2,1] * phix + tang[2,2]*phiy + normal[2] * pdou
    vel[3] = tang[3,1] * phix + tang[3,2]*phiy + normal[3] * pdou

    return vel
end

