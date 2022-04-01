using LinearAlgebra
using Interpolations
using DelimitedFiles
using Printf
using Plots

include("aeroel.jl")
include("geo.jl")
# See Katz Plotkin 12.1
# Lifting-Line Solution by horshoe elements

function test()

    max_iter = 500
    rlx = 0.4
    tol = 1e-6

    alpha = 2
    rho = 1.225
    V = 1

    # read in airfoil data
    c81 = readdlm("naca0012.csv", ',', Float64)
    cl_interp =  LinearInterpolation(c81[:,1], c81[:,2])
    cd_interp =  LinearInterpolation(c81[:,1], c81[:,3])

    span = 12
    chord = 1
    S = span.*chord
    npan = 8
    panVert, panCon, panCpt, bndLen, chordDir, panNorms = buildRectHShoeCCW(span, chord, npan)

    uinf = V.*[cosd(alpha), 0, sind(alpha)]

    useArtVisc = false

    A = zeros(Float64, npan, npan) # influence coefficents for normal flow
    B = zeros(Float64, npan, npan) # influence coefficents for downwash
    RHS = zeros(Float64, npan)
    alf = zeros(Float64, npan) # local angle of attack at each section
    w = zeros(Float64, npan) # downwash velocity
    ai = zeros(Float64, npan) # induced angle of attack
    alfe = zeros(Float64, npan) # effective angle of attack

    # preallocate working matrices
    X = zeros(2*npan) #(gam1, gam2, gamn ...., dalf1, dalf2, ....dalfn)
    F = zeros(2*npan) # nonlinear vector function
    J = zeros(2*npan, 2*npan)# jacobian of F

    # artificial viscosity
    mu = zeros(Float64, npan)
    #mu[:] .= 0.02

    # Influence matrix
    for i = 1:npan
        for j = 1:npan
            # influence of jth panel on ith collcation point
            uvw, dwash = hshoe(panVert[:,panCon[:,j]], panCpt[:,i])
            A[i,j] = dot(uvw, panNorms[:,i])
            B[i,j] = dot(dwash, panNorms[:,i])
        end
        RHS[i] = dot(-uinf, panNorms[:,i])
    end

    # calculate local alphas at eacch station
    for i = 1:npan
        vt = dot(uinf, chordDir[:,i])
        vp = dot(uinf, panNorms[:,i])
        alf[i] = atan(vp, vt)
    end

    # intial guess of X using angle of attack
    #cl_init = cl_interp(rad2deg.(alf))
    #gam_init = cl_init*chord*V/2

    gam_init = A\RHS

    println(gam_init)

    X[1:npan].= gam_init
    w[:] .= B*X[1:npan] #downash velocity
    ai[:] .= -atan.(w,V)# induced angle of attack

    done = false
    iter = 0

    # Find F(X)=0 using Newton Raphson

    while !done
        iter += 1
        for i = 1:npan
            # get neighboring panels
            if i == 1
                neigh = [i+1]
            elseif i ==npan
                neigh = [i-1]
            else
                neigh = [i-1; i+1]
            end

            # Lookup cl from angle of attack
            alfe[i] = rad2deg(alf[i]-ai[i]-X[npan+i])
            alfe[i] = deg180(alfe[i])
            clvisc = cl_interp(alfe[i])
            #clvisc = 2*pi*deg2rad(alfe[i]) # debug, assume 2pi airfoil

            # compute F(X)
            F[i] = sum(A[i,:].*X[1:npan]) + sin(alf[i]-X[npan+i])
            #F[i] = sum(A[i,:].*X[1:npan]) + sin(alf[i]-X[npan+i] - ai[i]) # ?
            F[npan+i] = X[npan+i]-(-2*X[i]/chord/V - clvisc)/(2*pi)

            # Compute J(X)
            for j = 1:npan
                J[i,j] = A[i,j]
            end
            J[i,npan+i] = -cos(alf[i]-X[npan+i])
            #J[i,npan+i] = -cos(alf[i]-X[npan+i]-ai[i]) # ?
            J[npan+i, i] = 2/chord/V /(2*pi)
            J[npan+i, npan+i] = 1

            # artificial viscosity contributions to F(X) and J(X)
            if useArtVisc
                dfdg = J[i,i]
                mu[i] = 1*(dfdg-1)/2
                mu[i] = 1e-2
                F[npan+i] = F[npan+i] +mu[i]*2*X[npan+i]
                J[npan+i,npan+i] = J[npan+i,npan+i] + 2*mu[i]
                for n = neigh
                    F[npan+i] = F[npan+i] -mu[i]*X[npan+n]
                    J[npan+i, npan+n] = J[npan+i, npan+n] - mu[i]
                end
            end
        end

        # calculate next iteration X
        Xnew = X - inv(J)*F.*rlx
        #println(maximum(abs.(F)))
        #println(Xnew[1:npan])

        # update X
        #println(rad2deg.(Xnew[npan+1:end]))
        X[:] .= Xnew
        w = B*X[1:npan] #downash velocity
        ai = -atan.(w,V)# induced angle of attack

        residual = maximum(abs.(F))
        if residual < tol
            done = true
        elseif iter == max_iter
            done = true
            println("Maximum iterations reached")
        end
    end
    
    # results
    res_cl = cl_interp(alfe)
    dL = 1/2 .*rho .* V.^2 .* bndLen .* chord .* res_cl
    dDi = -rho.*w.*X[1:npan].*bndLen
    L = sum(dL)
    Di = sum(dDi)
    CL = L/(1/2*rho*V^2*S)

    @printf "CL=%.8f\n" CL
    p1 = plot(panCpt[2,:], alfe)
    p2 = plot(panCpt[2,:], res_cl)
    p3 = plot(panCpt[2,:], 2 .*X[1:npan]./chord./V)
    display(p1)
    display(p2)
    #display(p3)
    println(alfe)
    println(rad2deg.(ai))
    println(rad2deg.(X[npan+1:end]))
    println(L)

    
end

test()