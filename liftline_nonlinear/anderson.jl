using Interpolations
using DelimitedFiles
using Plots

# Method from John D Anderson, Fundamentals of Aerodynamics

function test()

    c81 = readdlm("naca0012.csv", ',', Float64)
    cl_interp =  LinearInterpolation(c81[:,1], c81[:,2])
    cd_interp =  LinearInterpolation(c81[:,1], c81[:,3])

    a = 4;
    V = 10;
    n = 100
    v = 10 .*[cosd(a); 0; sind(a)] 
    chords = 1 .*ones(n+1)

    y = collect(LinRange(-10,10,n+1))
    b = y[end]-y[1]
    dy = diff(y)
    ymid = (y[2:end]+y[1:end-1])./2

    gam0 = 1
    gams = gam0.*sqrt.(1 .-(2 .*y./b).^2)
    new_gams = zeros(n+1)

    cls = zeros(n+1)
    alphas = zeros(n+1)

    done = false
    iter = 0

    while !done
        iter = iter + 1
        dgam = diff(gams)
        for i = 1:n+1
            ai = 0
            for j=1:n
                ai = ai + dgam[j]/dy[j]/(y[i]-ymid[j])*dy[j]
            end
            ai = 1/(4*pi*V)*ai
            alphas[i] = a - rad2deg(ai)
        end
        # extrapolate for alpha
        alphas[1] = (alphas[2]-alphas[3])/(y[2]-y[3])*(y[1]-y[3])+alphas[3]
        alphas[end] = (alphas[end-1]-alphas[end-2])/(y[end-1]-y[end-2])*(y[end]-y[end-2])+alphas[end-2]
        for i = 1:n+1
            cls[i] =  cl_interp(alphas[i])
            new_gams[i] = 1/2*V*chords[i]*cls[i]
        end
        # enforce circulation to be 0 at tips
        new_gams[1] = 0
        new_gams[end] = 0
       
        eps = abs.(new_gams-gams)
        if maximum(eps) < 1e-6
            done = true
        elseif iter == 500
            done = true
            println("Max iter reached")
        else
            k = 0.05
            gams = gams + k*(new_gams-gams)
        end

    end
    plot(y, cls)
end



test()
