using Interpolations
using DelimitedFiles
using Plots
using Printf

# strongly coupled algorithm from Gallay, Nonlinear Lifting Line for PrePost-stall Flows.
# Applied to a discrete wing modelled using vortex lattices
# Borrows geometry and visualization functions from VPM code

include("geo.jl")
include("vis.jl")
include("aeroel.jl")


function test()
    max_iter = 500
    rlx = 0.4
    tol = 1e-3

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
    npan = 80
    # read in panels
    panVert, panCon = createRect(span, chord, npan, 1,[0;0;0],0)
    panCpt, panNorm, panNPt, panEdgeVec, panEdgeLen, panEdgeUVec, panArea, panTang, panSinTi, panCosTi = calcPanProps(panVert, panCon)
    panNeighIdx, panNeighSide, panNeighDir, panNNeigh = calcneighbors(panCon, panNPt)


    npts = size(panVert,2)
    npan = size(panCon,2)

    panPres = zeros(Float32, npan)
    panVelSurf = zeros(Float32, 3, npan)
    panGam = zeros(Float64, npan)

    A = zeros(Float64, npan, npan)
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
        end
    end
    
    # build rhs vector
    for i = 1:npan
        RHS[i] = -dot(uinf, panNorm[:,i])
        #RHS[i] = RHS[i] - dot(panels[i].wake_vel, panNorm[:,i])
        #RHS[i] = RHS[i] + dot(panels[i].vcpt, panNorm[:,i])
    end
    
    # solve matrix for panel gamma
    #gam_init = A\RHS

    # calculate local alphas at each station
    for i = 1:npan
        alf[i] = deg2rad(alpha) 
    end

    # intial guess of X using angle of attack
    cl_init = cl_interp(rad2deg.(alf))
    gam_init = cl_init*chord*V/2

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
            F[npan+i] = X[npan+i]-(2*X[i]/chord/V - clvisc)/(2*pi)

            # Compute J(X)
            for j = 1:npan
                J[i,j] = A[i,j]
            end
            J[i,npan+i] = -cos(alf[i]-X[npan+i])
            J[npan+i, i] = -2/chord/V /(2*pi)
            J[npan+i, npan+i] = 1
        end

        display(J)
        quit()

        residual = maximum(abs.(F))
        if residual < tol
            done = true
        elseif iter == max_iter
            done = true
            println("Maximum iterations reached")
        end
    end

    pan2Vtu(panCon, panVert, panGam, panPres, panVelSurf, "test.vtu")
    quit()
end

test()