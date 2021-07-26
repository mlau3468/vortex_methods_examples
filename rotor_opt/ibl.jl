using DelimitedFiles
using LinearAlgebra


function panelMethod(x, y, s, theta, sine, cosine, rhs, xb, yb, alpha, bc)
    npts = length(xb)
    npan = npts - 1
    cn1 = similar(xb, npan, npan)
    cn2 = similar(xb, npan, npan)
    ct1 = similar(xb, npan, npan)
    ct2 = similar(xb, npan, npan)
    for i = 1:npan
        for j = 1:npan
            if i == j
                # kutta condition
                cn1[i,j] = -1
                cn2[i,j] = 1
                ct1[i,j] = pi/2
                ct2[i,j] = pi/2
            else
                a = -(x[i]-xb[j])*cosine[j] - (y[i]-yb[j])*sine[j]
                b = (x[i]-xb[j])^2 + (y[i]-yb[j])^2
                c = sin(theta[i]-theta[j])
                d = cos(theta[i]-theta[j])
                e = (x[i]-xb[j])*sine[j] - (y[i]-yb[j])*cosine[j]
                f = log(1+s[j]*(s[j]+2*a)/b)
                g = atan(e*s[j],b+a*s[j])
                p = (x[i]-xb[j])*sin(theta[i]-2*theta[j]) + (y[i]-yb[j])*cos(theta[i]-2*theta[j])        
                q = (x[i]-xb[j])*cos(theta[i]-2*theta[j]) - (y[i]-yb[j])*sin(theta[i]-2*theta[j])
                cn2[i,j] = d+0.5*q*f/s[j] - (a*c+d*e)*g/s[j]
                cn1[i,j] = 0.5*d*f +c*g -cn2[i,j]
                ct2[i,j] = c+0.5*p*f/s[j] + (a*d-c*e)*g/s[j]
                ct1[i,j] = 0.5*c*f-d*g-ct2[i,j]
            end
        end
    end
    # compute influence coefficients
    an = similar(xb, npan+1, npan+1)
    at = similar(xb, npan+1, npan+1)
    for i = 1:npan
        an[i,1] = cn1[i,1]
        an[i,npts] = cn2[i,npan]
        at[i,1] = ct1[i,1]
        at[i,npts] = ct2[i,npan]
        for j = 2:npan
            an[i,j] = cn1[i,j] + cn2[i,j-1]
            at[i,j] = ct1[i,j] + ct2[i, j-1]
        end
    end
    an[npts, 1] = 1
    an[npts, npts] = 1
    for j = 2:npan
        an[npts,j] = 0
    end
    rhs[npts] = 0
    gamma = inv(an)*(rhs+bc)

    v = similar(xb, npan)
    c_p = similar(xb, npan)
    for i = 1:npan
        v[i] = cos(theta[i] - alpha)
        for j = 1:npts
            v[i] = v[i] + at[i,j]*gamma[j]
            c_p[i] = (1-v[i]^2)
        end
    end
    return gamma, v, c_p
    
end

function airfoilCalc(pts, uinf, alf)
    alf = deg2rad(alf)
    xb = pts[:,1]
    yb = pts[:,2]
    npan = size(pts)[1]-1
    c = 1
    mu = 1.81206e-5
    rho = 1.225
    nu = mu/rho
    q = 1.2*rho*uinf^2
    re = uinf*c/nu
    threshold = 1e-5
    err = 1000
    iter = 0
    max_iter = 100

    # initialization
    x = similar(pts, npan)
    y = similar(pts, npan)
    s = similar(xb, npan)
    theta = similar(xb, npan)
    sine = similar(xb, npan)
    cosine = similar(xb, npan)
    rhs = similar(xb, npan+1)

    deltas = similar(pts, npan)
    thetas = similar(pts, npan)
    cf = similar(pts, npan)
    old_dels = similar(pts, npan)
    old_dels[:] .= 0
    g = similar(pts, npan+1)
    c_p = similar(pts, npan)
    vtan = similar(pts, npan)

    stag = 0
    trans = [0 0]
    sp = [0 0]

    g[:] .= 0.0

    # calculate panel properties
    for i = 1:npan
        ip1 = i+1
        x[i] = 0.5*(xb[i]+xb[ip1])
        y[i] = 0.5*(yb[i]+yb[ip1])
        s[i] = sqrt((xb[ip1]-xb[i])^2 + (yb[ip1]-yb[i])^2)
        theta[i] = atan(yb[ip1]-yb[i], xb[ip1]-xb[i])
        sine[i] = sin(theta[i])
        cosine[i] = cos(theta[i])
        rhs[i] = sin(theta[i]-alf)
    end

    
    while err >= threshold && iter < max_iter
        gamma, vtan, c_p = panelMethod(x, y, s, theta, sine, cosine, rhs, xb, yb, alf, g)
        ue = abs.(uinf*vtan)
        # get stagnation point
        stag = argmin(abs.(ue))
        up = stag:npan
        low = stag:-1:1
        # boundary layer solver
        deltas[up], thetas[up], cf[up], trans_u, sp_u = boundaryLayer(ue[up], x[up], y[up], c, nu)
        deltas[low], thetas[low], cf[low], trans_l, sp_l = boundaryLayer(ue[low], x[low], y[low], c, nu)

        # transition and separation points
        trans = [stag-1+trans_u stag+1-trans_l]
        sp = [stag-1+sp_u stag+1-sp_l]

        # iteration
        # calculate boundary condition for inviscid solver
        g[1:end-1] = 0.03*deltas
        g[end] = 0.03*deltas[npan]

        # looping criteria
        err = sum(abs.((deltas-old_dels)./old_dels)*100)
        old_dels[:] = deltas
        iter = iter + 1
        println(iter)
    end
    
    if iter < max_iter
        # trailing edge error handling
        c_p[1] = 1
        c_p[end] = 1
        vtan[1] = 0
        vtan[end] = 0
        # calculate cl, cd and plot
        e = 1.7
        cf = e*cf
        delta = e*deltas
        cp_u = 0
        cp_l = 0
        cd_lam = 0
        cd_turb = 0
        # upper airfoil
        for i = stag:npan-1
            cp_u = cp_u + (c_p[i+1] + c_p[i])*(x[i+1]-x[i])/2
            if i < trans[1]
                cd_lam = cd_lam + (cf[i+1] + cf[i])*(x[i+1]-x[i])/2
            else
                cd_turb = cd_turb + (cf[i+1] + cf[i])*(x[i+1]-x[i])/2
            end
        end
        # lower airfoil
            for i = stag:-1:2
                cp_l = cp_l + (c_p[i-1] + c_p[i])*(x[i-1]-x[i])/2
                if i > trans[2]
                    cd_lam = cd_lam + (cf[i-1] + cf[i])*(x[i-1]-x[i])/2
                else
                    cd_turb = cd_turb + (cf[i-1] + cf[i])*(x[i-1]-x[i])/2
                end
            end

        cl = (cp_l-cp_u)*cos(alf)
        cd = (cd_lam + cd_turb)
    else
        println("Convergence Failed.")
        cl = NaN
        cd = NaN
    end
    return cl, cd
end


function ft(t, cf, H, ue, due)
    return cf/2 - (2+H)*t*due/ue
end

function rungeKutta(theta0, x0, xStop, cf, H, ue, due)

    N = 2;
    theta = theta0
    thetap = theta0;
    
    eps = 1e-03
    temp = 1
    while temp > eps
        h = (xStop - x0)/N;
        for i = 1:N
            k1 = h*ft(theta, cf, H, ue, due)
            k1 = real(k1)
            k2 = h*ft(theta + k1/2, cf, H, ue, due)
            k2 = real(k2)
            k3 = h*ft(theta + k2/2, cf, H, ue, due)
            k3 = real(k3)
            k4 = h*ft(theta + k3, cf, H, ue, due)
            k4 = real(k4)
            theta = theta + (k1 + 2*k2 + 2*k3 + k4)/6
        end
        temp = abs(theta - thetap)
        thetap = theta
    end
    return theta
end

function fa(PI)
    kapa = 0.41
    return (2 + 3.179*PI + 1.5*PI^2)/(kapa*(1 + PI))
end

function re_t(PI, lambd, H)
    kapa = 0.41
    B = 5
    return (1 + PI)*exp(kapa*lambd - kapa*B - 2*PI)/(kapa*H)
end

function beta(PI)
    return -0.4 + 0.76*PI + 0.42*PI^2
end

function boundaryLayer(ue, xp, yp, c, nu)
    # initialization
    # momentum thickness calculation
    eps = 1e-5
    trans = NaN
    sp = NaN
    m = length(xp)
    thetas = similar(xp, m)
    # calcualte due/dx
    due = similar(xp, m)
    rex = 0
    for i =1:length(due)
        if i != length(ue)
            due[i] = (ue[i+1] - ue[i])/((xp[i+1] - xp[i])/c)
        else
            due[i] = (ue[i] - ue[i-1])/((xp[i] - xp[i-1])/c)
        end
    end
    
    # calcualte ds for integration
    rp = [xp yp]
    ds = similar(xp, m-1)
    for i =1:m-1
        ds[i] = abs(norm(rp[i+1,:] - rp[i,:]))
    end
    
    # laminar regime
    thetas[1] = sqrt(0.075*nu/abs(due[1]))
    # thetas for other points
    for i = 2:m
        integral = 0
        xt = 0
        for j = 2:i
            integral = integral + (ue[j]^5 + ue[j-1]^5)*ds[j-1]/2
            xt = xt + abs(xp[j] - xp[j-1])
        end
        thetas[i] = sqrt(0.45*nu*integral/(ue[i]^6))

        # checking transition point
        lhs = abs(ue[i])*thetas[i]/nu
        rex = abs(ue[i])*xt/nu
        # Cebeci and smith
        rhs = 1.174*(1+(22400/rex))*(rex^0.46)
        temp = abs(lhs-rhs)
        if isnan(trans*i)
            if temp < eps || lhs >= rhs
                trans = i
            end
        end
    end
    if isnan(trans) # if no transiiton found
        trans = m-1
    end
    # calculate lambda (pressure gradient parameter)
    lambda = (thetas.^2).*due/nu
    # calcualte tau_wall (wall shear stress) and deltas (disp. thickness)
    L = similar(xp, m)
    L[:] .= 0
    H = similar(xp, m)
    H[:] .= 0
    for i = 1:m
        z = 0.25 - lambda[i]
        if lambda[i] < 0.1 && lambda[i] > 0
            L[i] = 0.22 + 1.57*lambda[i] - 1.8*lambda[i]^2
            H[i] = 2.61 - 3.75*lambda[i] + 5.24*lambda[i]^2
        elseif lambda[i] <= 0 && lambda[i] > -0.1
            L[i] = 0.22 + 1.402*lambda[i] + (0.018*lambda[i]/(lambda[i]+0.107))
            H[i] = 2.088 + 0.0731/(lambda[i]+0.14)
        elseif lambda[i] >= 0.1 && lambda[i] <= 0.25
            if i == m
                L[i] = L[i-1]
            else
            L[i] = L[i+1]
            end
            H[i]= 2.0 + 4.14*z - 83.5*z^2 + 854*z^3 - 3337*z^4 + 4576*z^5
        else
            L[i] = L[i-1]
            H[i] = H[i-1]
        end
    end
    cf = 2*nu*L./(ue.*thetas)
    thetas = thetas
    cf[thetas .== 0] .= 0
    # turbulent regime
    e = 1
    for i = trans+1:m
        PI = 0.43
        cf0 = cf[i-1]
        temp = 1
        lam = sqrt(2/cf0)
        thetat = 0.0
        cft = 0.0
        while temp > eps
            Ht = lam/(lam-fa(PI))
            cft = 0.3*exp(-1.33*Ht)/(log10(re_t(PI,lam,Ht))^(1.74+0.31*Ht))

            thetat = rungeKutta(thetas[i-1], xp[i-1], xp[i], cft, Ht, ue[i], due[i])
            thetat = real(thetat)
            # iterator
            lam = sqrt(2/cft)
            betat = - Ht*thetat*due[i]*(lam^2)/ue[i]
            lam = lam + 0.1
            temp = abs(cft-cf0)
            cf0 = cft
        end
        test = isinf(thetat) || isnan(thetat)
        if test 
            thetas[i] = thetas[i-1]
            if isnan(i*sp)
                sp = i
            end
        else
            thetas[i] = e*thetat
        end
        cf[i] = e*cft
    end

    deltas = H.*thetas

    return real(deltas), real(thetas), real(cf), trans, sp
end
