using Plots: display
using Base: Float64
using Plots


function genNACA(n, nPts=100)
    nPts = nPts + 1
    c = 1
    ndig = length(n)
    # x coordinates
    cSpace = true
    if cSpace
        beta = LinRange(0, pi, nPts)
        x = (1 .-cos.(beta))/2
    else
        x = LinRange(0,1,nPts)
    end

    # thickness
    # maximum thickness as fraction of chord (last 2 digits)
    t = n[end-1]*10 + n[end]
    t = t/100
    openTE = false
    if openTE
        # thickness y coordinate with closed trailing edge
        y_t=t/0.2*(0.2969*sqrt.(x)-0.126*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4)
    else
        # thickness y coordinate with open trailing edge
        y_t=t/0.2*(0.2969*sqrt.(x)-0.126*x-0.3516*x.^2+0.2843*x.^3-0.1036*x.^4)
    end

    # camber
    y_c = zeros(nPts)
    dyc_dx = zeros(nPts)
    if ndig == 4 # NACA 4 series
        m = n[1]/100 # max camber
        p = n[2]/10 # lcoation of max camber
        if m == 0
            sym = 1
        else
            sym = 2
        end 

        for i=1:1:nPts
            if x[i] < p
                y_c[i]=m*x[i]/p^2*(2*p-x[i]) # mean camber y coordinate
                dyc_dx[i]=2*m/p^2*(p-x[i]) # mean camber first derivative
            else
                y_c[i]=m*(1-x[i])/(1-p)^2*(1+x[i]-2*p)
                dyc_dx[i]=2*m/(1-p)^2*(p-x[i])
            end
        end
    elseif  ndig == 5 # NACA 5 series
        p = n[2]/20 #location of maximum camber
        rn = n[3] # type of camber
        if rn == 0 # standard camber
            r=3.33333333333212*p^3+0.700000000000909*p^2+1.19666666666638*p-0.00399999999996247
            k1=1514933.33335235*p^4-1087744.00001147*p^3+286455.266669048*p^2-32968.4700001967*p+1420.18500000524
            k2_k1=85.5279999999984*p^3-34.9828000000004*p^2+4.80324000000028*p-0.21526000000003
            for i = 1:1:nPts
                if x[i]< r
                    y_c[i]=k1/6*(x[i]^3-3*r*x[i]^2+r^2*(3-r)*x[i])
                    dyc_dx[i]=k1/6*(3*x[i]^2-6*r*x[i]+r^2*(3-r))
                else
                    y_c[i]=k1*r^3/6*(1-x[i])
                    dyc_dx[i]=-k1*r^3/6
                end
            end
        elseif  rn == 1 # reflexed camber
            r=10.6666666666861*p^3-2.00000000001601*p^2+1.73333333333684*p-0.0340000000002413
            k1=-27973.3333333385*p^3+17972.8000000027*p^2-3888.40666666711*p+289.076000000022
            k2_k1=85.5279999999984*p^3-34.9828000000004*p^2+4.80324000000028*p-0.21526000000003
        
            for i = 1:1:nPts
                if x[i]< r
                    y_c[i]=k1/6*((x[i]-r)^3-k2_k1*(1-r)^3*x[i]-r^3*x[i]+r^3)
                    dyc_dx[i]=k1/6*(3*(x[i]-r)^2-k2_k1*(1-r)^3-r^3)
                else
                    y_c[i]=k1/6*(k2_k1*(x[i]-r)^3-k2_k1*(1-r)^3*x[i]-r^3*x[i]+r^3)
                    dyc_dx[i]=k1/6*(3*k2_k1*(x[i]-r)^2-k2_k1*(1-r)^3-r^3)
                end
            end
        end
    elseif ndig == 6 # 6 digit
        ser = n[1] # number of series
        a = n[2]/10 # chordwise position of minimum pressure
        c_li = n[4]/10 # design lift coefficient
        g=-1/(1-a)*(a^2*(1/2*log(a)-1/4)+1/4)
        h=1/(1-a)*(1/2*(1-a)^2*log(1-a)-1/4*(1-a)^2)+g
        if ser == 6
            y_c=c_li/(2 .*pi.*(a.+1)).*(1/(1 .-a).*(1/2 .*(a .-x).^2 .*log.(abs.(a .-x))-1/2 .*(1 .-x).^2 .*log.(1 .-x) .+1/4 .*(1 .-x).^2-1/4 .*(a .-x).^2)-x.*log.(x) .+g .-h.*x)
            dyc_dx=.-(c_li.*(h.+log.(x)-(x/2 .-a/2+(log.(1 .-x).*(2 .*x.-2))/2 .+(log.(abs.(a .-x)).*(2*a.-2*x))/2+(sign.(a.-x).*(a.-x).^2)./(2*abs.(a.-x)))/(a.-1).+1))/(2*pi*(a .+1))
        end

    end 

    # final calculations
    theta = atan.(dyc_dx) # angle for modifying x coordinate
    # assign coordinates
    x_e = (x-y_t.*sin.(theta))*c
    x_i = (x+y_t.*sin.(theta))*c
    y_e = (y_c+y_t.*cos.(theta))*c
    y_i = (y_c-y_t.*cos.(theta))*c
    if ndig == 6
        x_e[1] = 0
        x_e[end] = 1
        x_i[1] = 0
        x_i[end] = 1
        y_e[1] = 0
        y_e[end] = 0
        y_i[1] = 0
        y_i[end] = 0
    end


    # to selig format
    x_i = x_i[end:-1:1]
    y_i = y_i[end:-1:1]
    
    xcoord = vcat(x_i, x_e[2:end])
    ycoord = vcat(y_i, y_e[2:end])

    pts = [xcoord ycoord]

    return pts
end

