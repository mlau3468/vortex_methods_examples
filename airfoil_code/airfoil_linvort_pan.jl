include("airfoil_util.jl")

function linVortPan(mu1, mu2, p1, p2, p, coloc=false)
    r1 = ptDist(p1, p)
    r2 = ptDist(p2, p)
    theta = -atan(p2[2]-p1[2], p2[1]-p1[1])
    p_new = coordRot2D(p, theta, p1)

    x1 = 0
    x2 = ptDist(p1, p2)

    x = p_new[1]
    z = p_new[2]

    theta1 = atan(z,x)
    theta2 = atan(z, x-x2)

    tol = 1e-8
    # if on panel itself
    if (z < tol && abs(theta1) < tol && abs((abs(theta2)-pi) < tol)) || coloc
        up_a=-0.5*(x-x2)/(x2)
        up_b=0.5*(x)/(x2)
        wp_a=-1/2/pi
        wp_b=1/2/pi
    else
        # first panel point
        up_a = z/2/pi*(-mu1)/(x2-x1)*log(r2/r1) + (mu1*(x2-x1)-mu1*(x-x1)) /2/pi/(x2-x1) * (theta2-theta1)
        wp_a = -(mu1*(x2-x1)-mu1*(x-x1)) /2/pi/(x2-x1) * log(r1/r2) + z/2/pi*(-mu1)/(x2-x1) * ((x2-x1)/z-(theta2-theta1))
        # katz plotkin pg 303. note error in wp term. a + is actually a -
        # second panel point
        up_b = z/2/pi*(mu2)/(x2-x1)*log(r2/r1) + mu2*(x-x1)/2/pi/(x2-x1) * (theta2-theta1)
        wp_b = -(mu2*(x-x1))/2/pi/(x2-x1)*log(r1/r2) + z/2/pi*(mu2)/(x2-x1) * ((x2-x1)/z-(theta2-theta1))
    end
    uw_a = coordRot2D([up_a, wp_a], -theta, [0,0])
    uw_b = coordRot2D([up_b, wp_b], -theta, [0,0])

    uw = uw_a .+ uw_b
    
    return [uw'; uw_a'; uw_b']

end

function airfoil_linVort(pan_pts, alpha)
pan_pts, c_pts, thetas, norms, tangents, dists = procPanels(pan_pts)
# conditions
U = 1
chord = 1
rho = 1.225

# Initialize solver matrix
num_pan = size(pan_pts,1)
A = zeros(num_pan+1, num_pan+1) # normal matrix
B = zeros(num_pan, num_pan+1) # tangent matrix
RHS = zeros(num_pan+1)
u_vec = U .* [cosd(alpha), sind(alpha)]

for i = 1:num_pan
    RHS[i] = -u_vec'norms[i, :]
    for j = 1:num_pan
        if i == j 
            res = linVortPan(1, 1, pan_pts[j,1:2], pan_pts[j,3:4], c_pts[i,:], true)
        else
            res = linVortPan(1, 1, pan_pts[j,1:2], pan_pts[j,3:4], c_pts[i,:])
        end 
        
        uw_a = res[2,:]
        uw_b = res[3,:]
        
        if j == 1
            A[i,j] = uw_a'norms[i,:]
            B[i,j] = uw_b'tangents[i,:]
        elseif j == num_pan
            A[i,j] = uw_a'norms[i,:] + holda
            B[i,j] = uw_b'tangents[i,:] + holdb
            A[i,num_pan+1] = uw_b'norms[i,:]
            B[i,num_pan+1] = uw_b'tangents[i,:]
        else
            A[i,j] = uw_a'norms[i,:] + holda
            B[i,j] = uw_b'tangents[i,:] + holdb
        end

        global holda = uw_b'norms[i,:]
        global holdb = uw_b'tangents[i,:]
    end
end

# Kutta condition
A[end, 1] = 1
A[end, end] = 1
RHS[end] = 0

mu = A\RHS

# lift by kutta jokowski
ls = zeros(num_pan)
for i = 1:num_pan
    ls[i] = rho*U.*(mu[i] + mu[i+1])/2*dists[i]
end
lift = sum(ls)
cl = lift/0.5/rho/U^2/chord

# get cp
vs = zeros(num_pan)
cp = zeros(num_pan)
for i = 1:num_pan
    vs[i] = B[i,:]'mu + u_vec'tangents[i,:]
    cp[i] = 1-vs[i]^2/U^2
end
cpPlot = plot(c_pts[:,1], cp, yflip=true)

return [cl, cpPlot]
end



# Testing
#pan_pts = readAF("airfoil.csv", true)
pan_pts = readAF("4521.csv", true)
#pan_pts = repanel(pan_pts, 80, 0.75, true)

cl, cpPlot = airfoil_linVort(pan_pts, 2)
display(cpPlot)
println("CL: $cl")