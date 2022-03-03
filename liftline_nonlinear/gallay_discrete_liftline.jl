using LinearAlgebra
using Interpolations
using DelimitedFiles
using Printf
using Plots

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
    chordDir = zeros(3,n)

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


function test()
    # read in airfoil data
    c81 = readdlm("naca0012.csv", ',', Float64)
    cl_interp =  LinearInterpolation(c81[:,1], c81[:,2])
    cd_interp =  LinearInterpolation(c81[:,1], c81[:,3])

    span = 12
    chord = 1
    S = span.*chord
    npan = 10
    panVerts, panCon, panCpts, bndLen, chordDir = buildRectHShoe(span, chord, npan)
    panNorms = calcPanNorm(panVerts, panCon)

    alpha = 5
    rho = 1.225
    V = 1
    uinf = V.*[cosd(alpha), 0, sind(alpha)]

    A = zeros(Float64, npan, npan) # influence coefficents for normal flow
    B = zeros(Float64, npan, npan) # influence coefficents for downwash
    RHS = zeros(Float64, npan)
    alf = zeros(Float64, npan) # local angle of attack at each section
    w = zeros(Float64, npan) # downwash velocity
    ai = zeros(Float64, npan) # induced angle of attack

    # preallocate working matrices
    X = zeros(2*npan) #(gam1, gam2, gamn ...., dalf1, dalf2, ....dalfn)
    F = zeros(2*npan) # nonlinear vector function
    J = zeros(2*npan, 2*npan)# jacobian of F


    for i = 1:npan
        for j = 1:npan
            # influence of jth panel on ith collcation point
            uvw, dwash = hshoe(panVerts[:,panCon[1,j]], panVerts[:,panCon[2,j]], panVerts[:,panCon[3,j]], panVerts[:,panCon[4,j]], panCpts[:,i])
            A[i,j] = dot(uvw, panNorms[:,i])
            B[i,j] = dot(dwash, panNorms[:,i])
        end
        #RHS[i] = -dot(uinf, panNorms[:,i])
        RHS[i] = -sind(alpha)
    end

    # intial guess of X using uncorrected gamma
    X[1:npan].= A\RHS
    w[:] .= B*X[1:npan] #downash velocity
    ai[:] .= -atan.(w,V)# induced angle of attack

    # calculate local alphas at eacch station
    for i = 1:npan
        alf[i] = deg2rad(alpha) 
    end

    done = false
    max_iter = 100
    iter = 0
    tol = 1e-6

    # Find F(X)=0 using Newton Raphson

    while !done
        iter += 1
        # compute F(X)
        for i = 1:npan
            #F[i] = sum(A[i,:].*X[1:npan]) - sin(alf[i]-X[npan+i])
            F[i] = sum(A[i,:].*X[1:npan]) + sin(alf[i]-X[npan+i])
            
            alfe = rad2deg(alf[i]-ai[i]-X[npan+i])
            alfe = mod(alfe, 360)
            if alfe > 180
                alfe = alfe - 360
            elseif alfe < -180
                alfe = alfe + 360
            end

            clvisc = cl_interp(alfe)
            #F[npan+i] = X[npan+i]-(clvisc - 2*X[i]/chord/V)/(2*pi)
            F[npan+i] = X[npan+i]-(2*X[i]/chord/V - clvisc)/(2*pi)
        end

        # compute J(X)
        for i =1:npan
            for j = 1:npan
                J[i,j] = A[i,j]
            end
            #J[i,i+npan] = cos(alf[i]-X[npan+i])
            J[i,i+npan] = -cos(alf[i]-X[npan+i])
            #J[i+npan, i] = 2/chord/V /(2*pi)
            J[i+npan, i] = -2/chord/V /(2*pi)
            J[i+npan, i+npan] = 1
        end

        # calculate next iteration X
        Xnew = X - inv(J)*F

        residual = maximum(abs.(Xnew-X))
        if residual < tol
            done = true
        elseif iter == max_iter
            done = true
            println("Maximum iterations reached")
        end

        # update X
        #println(rad2deg.(Xnew[npan+1:end]))
        X[:] .= Xnew
        w = B*X[1:npan] #downash velocity
        ai = -atan.(w,V)# induced angle of attack
    end
    
    dL = rho.*V.*X[1:npan].*bndLen
    dDi = -rho.*w.*X[1:npan].*bndLen
    L = sum(dL)
    Di = sum(dDi)
    CL = L/(1/2*rho*V^2*S)

    @printf "CL=%.8f" CL
    plot(panCpts[2,:], X[1:npan])

    
end

test()