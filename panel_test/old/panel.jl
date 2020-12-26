using HDF5
using Statistics
using LinearAlgebra
using DelimitedFiles
using WriteVTK

struct refFrame
    name :: String
    parent :: String
    children # list of child references
    move2Parent :: Bool # Is this reference frame moving with respect to the parent frame
    move2Global :: Bool # Is this reference frame moving with respect to the global frame
    origin # 3x1 vector origin in the parent frame
    orient # 3x3 orientation matrix with respect to parent

    vel_g # 3x1 vector velocity with respect to global frame 
    rot_g #3x3 rotation rate matrix with respect to global frame
end

struct component
    name::String
    elType :: String
    refId :: Int
    refName :: String
    ee  # panels
    rr  # vertices
    neigh # neighboring
    n_pan :: Int # Number of panels
    n_vert :: Int # Number of vertices

    e_te #panel index associated with each TE panel
    i_te #component panel id of the nodes on the TE
    ii_te #TE id of the nodes of the TE elements (index in rr_te)
    rr_te #Coordinates of the nodes of the TE
    neigh_te #TE id of neighboring TE elements
    n_pan_te :: Int
    n_vert_te :: Int
    dir_te #Unit vector at TE nodes

end 

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

function sourc(pts, loc, dou, normal)
    
    cpt = mean(pts, dims=2)
    radius = norm(loc.-cpt)
    e3 = normal
    zQ = sum((loc.-cpt).*e3)
    Qp = loc .- zQ .* e3

    # get edge lens
    edge_len = zeros(4)
    edge_len[1] = norm(pts[:,2] .- pts[:,1])
    edge_len[2] = norm(pts[:,3] .- pts[:,2])
    edge_len[3] = norm(pts[:,4] .- pts[:,3])
    edge_len[4] = norm(pts[:,1] .- pts[:,4])
    # get edge vectors
    edge_vec = zeros(3,4)
    edge_vec[:,1] = pts[:,2] .- pts[:,1]
    edge_vec[:,2] = pts[:,3] .- pts[:,2]
    edge_vec[:,3] = pts[:,4] .- pts[:,3]
    edge_vec[:,4] = pts[:,1] .- pts[:,4]

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
    n_ver = 4
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


function dub(pts, loc, normal)

    cpt = mean(pts, dims=2)
    radius = norm(loc.-cpt)
    e3 = normal

    # get edge vectors
    edge_vec = zeros(3,4)
    edge_vec[:,1] = pts[:,2] .- pts[:,1]
    edge_vec[:,2] = pts[:,3] .- pts[:,2]
    edge_vec[:,3] = pts[:,4] .- pts[:,3]
    edge_vec[:,4] = pts[:,1] .- pts[:,4]

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


function eerr2vtk(ee, rr, fname)
    celltype = VTKCellTypes.VTK_QUAD
    cells = MeshCell[]
    cdata = Float32[]
    pts = rr
    for i =1:size(ee,2)
        inds = ee[:,i]
        c = MeshCell(celltype, inds)
        push!(cells, c)
    end
    outfile = vtk_grid(fname, pts, cells, compress=2) do vtk
    end
end

function vel_dub(pts, loc, normal)
    #=
    !> compute velocity AIC of panel <this>, on the control point <pos>
!! Biot-Savart law. Regular kernel: linear core (Rankine vortex)
!! Compute the velocity induced by a vortex ring (equivalent to a constant
!!  intentsity surface doublet) with intensity 4*pi. <------
=#
    v_dou = 0.0
    n_ver = 4

    # get edge vectors
    edge_vec = zeros(3,4)
    edge_vec[:,1] = pts[:,2] .- pts[:,1]
    edge_vec[:,2] = pts[:,3] .- pts[:,2]
    edge_vec[:,3] = pts[:,4] .- pts[:,3]
    edge_vec[:,4] = pts[:,1] .- pts[:,4]
    # get edge lens
    edge_len = zeros(4)
    edge_len[1] = norm(pts[:,2] .- pts[:,1])
    edge_len[2] = norm(pts[:,3] .- pts[:,2])
    edge_len[3] = norm(pts[:,4] .- pts[:,3])
    edge_len[4] = norm(pts[:,1] .- pts[:,4])
    # unit vectors
    edge_uni = zeros(3,4)
    edge_uni[:,1] = edge_vec[:,1] ./ norm(edge_vec[:,1])
    edge_uni[:,2] = edge_vec[:,2] ./ norm(edge_vec[:,2])
    edge_uni[:,3] = edge_vec[:,3] ./ norm(edge_vec[:,3])
    edge_uni[:,4] = edge_vec[:,4] ./ norm(edge_vec[:,4])

    r_rankine = 0.1
    r_cutoff = 0.001
    for i = 1:n_ver
        indp1 = 1 + i % n_ver
        indm1 = n_ver - (n_ver-i+1) % n_ver
        av = pos-pts[:i]
        ai = sum(av.*edge_uni(:,i))
        R1 = norm(av)
        R2 = norm(loc./pts[:,indp1])
        hv = av .- ai.*edge_uni[:,i]
        hi = norm(hv)

        if hi > edge_len[i]*r_rankine
            global v_dou = v_dou + ((edge_len[i]-ai)/r2 + ai/r1) / (hi*2) * cross(edge_uni[:,i], hv)
        else
            if R1 > edge_len[i]*r_cutoff && R2 > edge_len[i]*r_cutoff
                r_Ran = r_rankine * edge_len[i]
                global v_dou = v_dou + ((edge_len[i]-ai)/r2 + ai/r1) / (r_Ran*2) * cross(edge_uni[:,i], hv)
            end
        end
    end
    return v_dou
end

function vel_sourc(pts, loc, normal)
    phix = 0.0
    phiy = 0.0

    n_ver = 4

    # settings
    prev_qua = [4 1 2 3]
    next_qua = [2 3 4 1]
    prev_tri = [3 1 2]
    next_tri = [2 3 1]

    # get edge vectors
    edge_vec = zeros(3,4)
    edge_vec[:,1] = pts[:,2] .- pts[:,1]
    edge_vec[:,2] = pts[:,3] .- pts[:,2]
    edge_vec[:,3] = pts[:,4] .- pts[:,3]
    edge_vec[:,4] = pts[:,1] .- pts[:,4]
    # get edge lens
    edge_len = zeros(4)
    edge_len[1] = norm(pts[:,2] .- pts[:,1])
    edge_len[2] = norm(pts[:,3] .- pts[:,2])
    edge_len[3] = norm(pts[:,4] .- pts[:,3])
    edge_len[4] = norm(pts[:,1] .- pts[:,4])

    # unit vectors
    edge_uni = zeros(3,4)
    edge_uni[:,1] = edge_vec[:,1] ./ norm(edge_vec[:,1])
    edge_uni[:,2] = edge_vec[:,2] ./ norm(edge_vec[:,2])
    edge_uni[:,3] = edge_vec[:,3] ./ norm(edge_vec[:,3])
    edge_uni[:,4] = edge_vec[:,4] ./ norm(edge_vec[:,4])

    # central point
    cpt = mean(pts, dims=2)

    # local tangent unit vector as in PANAIR
    n_sides = 4
    tanl = 0.5 .* (pts[:,n_sides] .+ pts[:,1]) .- cpt
    tang = zeros(3,2)
    tang[:,1] = tanl ./ norm(tanl)
    tang[:,2] = cross(normal, tang[:,1])

    sinTi = zeros(n_sides)
    cosTi = zeros(n_sides)

    for i = 1:n_sides
        sinTi[i] = sum(edge_uni[:,i] .* tang(:,1))
        cosTi[i] = sum(edge_uni[:,i] .* tang(:,2))
    end


    for i = 1:n_ver
        if n_ver == 3 # triangle element
            indm1 = prev_tri[i]
            indp1 = next_tri[i]
        elseif n_ver == 4 # quad element
            indm1 = prev_qua[i]
            indp1 = next_qua[i]
        end

        R1 = norm(loc .- pts(:,i))
        R2 = norm(loc .- pts(:, indp1))

        if (R1+R2-edge_len[i < 1e-12])
            souLog = 0
        else
            souLog = log((R1+R2+edge_len[i]) / (R1+R2-edge_len[i]))
        end
        phix = phix + sinTi[i] * souLog
        phiy = phiy - cosTi[i] * souLog
    end

    pdou = dub(pts, loc, normal)

    phix = - phix
    phiy = - phiy
    pdou = - pdou

    vel = zeros(3)
    vel[1] = tang[1,1] * phix + tang[1,2]*phiy + normal[1] * pdou
    vel[2] = tang[2,1] * phix + tang[2,2]*phiy + normal[2] * pdou
    vel[3] = tang[3,1] * phix + tang[3,2]*phiy + normal[3] * pdou

    return vel
end


function hello(pts, loc, normal)

    cpt = mean(pts, dims=2)
    radius = norm(loc.-cpt)
    e3 = normal

    # get edge vectors
    edge_vec = zeros(3,4)
    edge_vec[:,1] = pts[:,2] .- pts[:,1]
    edge_vec[:,2] = pts[:,3] .- pts[:,2]
    edge_vec[:,3] = pts[:,4] .- pts[:,3]
    edge_vec[:,4] = pts[:,1] .- pts[:,4]
    
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