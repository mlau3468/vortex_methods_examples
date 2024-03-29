using LinearAlgebra
# See Katz Plotkin 12.1
# Lifting-Line Solution by horshoe elements

function vortxl!(x,y,z,x1,y1,z1,x2,y2,z2, gama, vel)
    #CALCULATES THE INDUCED VELOCITY (U,V,W) AT A POI
    #(X,Y,Z) DUE TO A VORTEX ELEMENT VITH STRENGTH GAMA PER UNIT LENGTH
    #POINTING TO THE DIRECTION (X2,Y2,Z2)-(X1,Y1,Z1).
    # r1x21
    rcut = 1e-10
    r1r2x = (y-y1)*(z-z2)-(z-z1)*(y-y2)
    r1r2y = -((x-x1)*(z-z2)-(z-z1)*(x-x2))
    r1r2z = (x-x1)*(y-y2)-(y-y1)*(x-x2)
    # (r1xr2)^2
    square = r1r2x*r1r2x+r1r2y*r1r2y+r1r2z*r1r2z
    # R0(R1/R(R1)-R2/R(R2))
    r1=sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1))
    r2=sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2))
    if (r1<rcut) || (r2<rcut) || (square<rcut)
        #vel[1] += 0
        #vel[2] += 0
        #vel[3] += 0
    else
        r0r1=(x2-x1)*(x-x1)+(y2-y1)*(y-y1)+(z2-z1)*(z-z1)
        r0r2=(x2-x1)*(x-x2)+(y2-y1)*(y-y2)+(z2-z1)*(z-z2)
        coef=gama/(4.0*pi*square)*(r0r1/r1-r0r2/r2)
        vel[1] += r1r2x*coef
        vel[2] += r1r2y*coef
        vel[3] += r1r2z*coef
    end

end

function hshoe(p1, p2, p3, p4, p)
    # calculates influence from one horshoe vortex on point p
    # pts A,B,C,D
    # vortexlines: A-B, B-C, C-D 
    uvw = zeros(3)
    vortxl!(p[1],p[2],p[3],p1[1],p1[2],p1[3],p2[1],p2[2],p2[3], 1.0, uvw)
    vortxl!(p[1],p[2],p[3],p2[1],p2[2],p2[3],p3[1],p3[2],p3[3], 1.0, uvw)
    vortxl!(p[1],p[2],p[3],p3[1],p3[2],p3[3],p4[1],p4[2],p4[3], 1.0, uvw)

    dwash = zeros(3)
    vortxl!(p[1],p[2],p[3],p1[1],p1[2],p1[3],p2[1],p2[2],p2[3], 1.0, dwash)
    vortxl!(p[1],p[2],p[3],p3[1],p3[2],p3[3],p4[1],p4[2],p4[3], 1.0, dwash)
    return uvw, dwash
end

function buildRectHShoe(span, chord, n)
    dy = span/n

    panVerts = zeros(Float64, 3, 2*(n+1))
    panCon = zeros(Int64,4, n)
    panCpts = zeros(Float64, 3, n)# collocation points
    bndLen = zeros(n)

    le_loc = 0.25

    npt = 1
    for i = 1:n
        p1 = [span*20; dy*(i-1); 0]
        p2 = [le_loc*chord; dy*(i-1); 0]
        p3 = [le_loc*chord; dy*(i); 0]
        p4 = [span*20; dy*(i); 0]

        # panel te for collcation point calculation
        p1c = [(le_loc+1)*chord; dy*(i-1); 0]
        p4c = [(le_loc+1); dy*(i); 0]

        # bound vortex segment length
        bndLen[i] = norm(p3-p2)

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

    return panVerts, panCon, panCpts, bndLen
end

function calcPanNorm(panVerts, panCon)
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


function run()
    span = 12
    chord = 1
    S = span.*chord
    npan = 10
    panVerts, panCon, panCpts, bndLen = buildRectHShoe(span, chord, npan)
    panNorms = calcPanNorm(panVerts, panCon)

    alpha = 5
    rho = 1.225
    V = 1
    uinf = V.*[cosd(alpha), 0, sind(alpha)]

    A = zeros(Float64, npan, npan) # influence coefficents for normal flow
    B = zeros(Float64, npan, npan) # influence coefficents for downwash
    RHS = zeros(Float64, npan)
    for i = 1:npan
        for j = 1:npan
            # influence of jth panel on ith collcation point
            uvw, dwash = hshoe(panVerts[:,panCon[1,j]], panVerts[:,panCon[2,j]], panVerts[:,panCon[3,j]], panVerts[:,panCon[4,j]], panCpts[:,i])
            A[i,j] = dot(uvw, panNorms[:,i])
            B[i,j] = dot(dwash, panNorms[:,i])
        end
        RHS[i] = -dot(uinf, panNorms[:,i])
    end

    gam = A\RHS
    display(gam)
    w = B*gam #downash velocity
    display(w)

    dL = rho.*V.*gam.*bndLen
    dDi = -rho.*w.*gam.*bndLen

    L = sum(dL)
    Di = sum(dDi)

    CL = L/(1/2*rho*V^2*S)

    
end

run()