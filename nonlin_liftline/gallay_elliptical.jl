using Interpolations
using DelimitedFiles
using Plots
using Printf

# strongly coupled algorithm from Gallay, Nonlinear Lifting Line for PrePost-stall Flows

function clvi(alpha)
    return pi*sin(2*alpha)
end

function clinv(alpha, da, lamda)
    # da = alpha correction factor
    # lambda = aspect ratio
    return 2*pi*(alpha-da)/(1+2/lamda)
end

function test()
    AR = 12.75
    alf = 20 # degrees

    tol = 1e-10

    A = zeros(2,2)
    B = zeros(2,1)
    sol = zeros(2,1) #(ae, da)
    alf = deg2rad(alf)

    iter = 0

    done = false
    while !done
        iter = iter + 1
        ae = sol[1]
        da = sol[2]
        A[1,1] = 1
        A[1,2] = 1/(1+2/AR)-1
        A[2,1] = 1-cos(2*ae)
        A[2,2] = -1

        B[1,1] = alf/(1+2/AR)
        B[2,1] = 1/2*sin(2*ae)-ae*cos(2*ae)

        x = A\B
        if maximum(abs.(sol-x)) < tol
            done = true
        elseif iter == 100
            done = true
            println("Max iter reached")
        end
        sol[:] .= x[:]
        println(sol)
    end
    # calculate cl
    cl = clvi(sol[1])
    @printf "CL=%.4f" cl

end


test()