using Interpolations
using DelimitedFiles
using Plots
using Printf
#=
strongly coupled algorithm from Gallay, Nonlinear Lifting Line for PrePost-stall Flows.
Applied to a discrete wing modelled using vortex rings. The wake is not modelled, so this will not 
yield a correct solution. Only meant to demonstrate one convergence loop for the 1st timestep in an
unsteady simulation.
Borrows geometry and visualization functions from VPM code
=#

include("geo.jl")
include("vis.jl")
include("aeroel.jl")


function test()
    max_iter = 500
    rlx = 0.4
    tol = 1e-6

    alpha = 2
    rho = 1.225
    V = 1

    useArtVisc = false
    uinf = V.*[cosd(alpha), 0, sind(alpha)]

    # read in airfoil data
    c81 = readdlm("naca0012.csv", ',', Float64)
    cl_interp =  LinearInterpolation(c81[:,1], c81[:,2])
    cd_interp =  LinearInterpolation(c81[:,1], c81[:,3])

    span = 12
    chord = 1
    S = span.*chord
    npan = 8
    # read in panels
    panVert, panCon, panCpt, bndLen, chordDir, panNorm, panArea = buildRect(span, chord, npan)

    panPres = zeros(Float32, npan)
    panVelSurf = zeros(Float32, 3, npan)
    panGam = zeros(Float64, npan)

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

    # Influence matrix
    for i = 1:npan
        for j = 1:npan
            # influence of jth panel on ith collocation point
            vel = vrtxring(panVert[:,panCon[:,j]], panCpt[:,i], 1.0)
            A[i,j] = dot(vel, panNorm[:,i])
            dwash = vrtxline(panVert[:,panCon[1,j]], panVert[:,panCon[2,j]], panCpt[:,i], 1)
            dwash = dwash + vrtxline(panVert[:,panCon[3,j]], panVert[:,panCon[4,j]], panCpt[:,i], 1)
            B[i,j] = dot(dwash, panNorm[:,i])
        end
        RHS[i] = dot(-uinf, panNorm[:,i])
    end

    # calculate local alphas at each station
    for i = 1:npan
        vt = dot(uinf, chordDir[:,i])
        vp = dot(uinf, panNorm[:,i])
        alf[i] = atan(vp, vt)
    end

    # solve matrix for panel gamma
    gam_init = A\RHS

    # intial guess of X using angle of attack
    #cl_init = cl_interp(rad2deg.(alf))
    #gam_init = cl_init*chord*V/2

    println(gam_init)

    done = false
    iter = 0

    # Find F(X)=0 using Newton Raphson

    while !done
        iter += 1
        for i = 1:npan
            # Lookup cl from angle of attack
            alfe[i] = rad2deg(alf[i]-ai[i]-X[npan+i])
            alfe[i] = deg180(alfe[i])
            clvisc = cl_interp(alfe[i])

            # compute F(X)
            F[i] = sum(A[i,:].*X[1:npan]) + sin(alf[i]-X[npan+i])
            F[npan+i] = X[npan+i]-(-2*X[i]/chord/V - clvisc)/(2*pi)

            # Compute J(X)
            for j = 1:npan
                J[i,j] = A[i,j]
            end
            J[i,npan+i] = -cos(alf[i]-X[npan+i])
            J[npan+i, i] = 2/chord/V /(2*pi)
            J[npan+i, npan+i] = 1
        end

        # calculate next iteration X
        Xnew = X - inv(J)*F.*rlx
        #println(maximum(abs.(F)))
        #println(Xnew[1:npan])

        # update X
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

    println(alfe)
    println(rad2deg.(ai))
    println(rad2deg.(X[npan+1:end]))

    # results
    res_cl = cl_interp(alfe)
    dL = 1/2 .*rho .* V.^2 .* panArea .* res_cl
    #dDi = -rho.*w.*X[1:npan].*bndLen
    L = sum(dL)
    #Di = sum(dDi)
    CL = L/(1/2*rho*V.^2*S)
    
    @printf "CL=%.8f\n" CL
    p1 = plot(panCpt[2,:], alfe)
    p2 = plot(panCpt[2,:], res_cl)
    p3 = plot(panCpt[2,:], -2 .*X[1:npan]./chord./V)
    display(p1)
    display(p2)
    display(p3)

    pan2Vtu(panCon, panVert, panGam, panPres, panVelSurf, "test.vtu")

end

test()