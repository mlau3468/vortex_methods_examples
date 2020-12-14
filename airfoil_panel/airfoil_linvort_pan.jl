include("airfoil_util.jl")

function linVort(mu1, mu2, p1, p2, p)
    r1 = ptDist(p1, p)
    r2 = ptDist(p2, p)
    theta = -atan(p2[2]-p1[2], p2[1]-p1[1]);
    p_new = coordRot2D(p, theta, p1)

    x1 = 0
    x2 = ptDist(p1, p2)

    x = p_new[1]
    z = p_new[2]
    theta1 = atan(z,x)
    theta2 = atan(z, x-x2)

    # first panel point
    up_a = z/2/pi*(-mu1)/(x2-x1)*log(r2/r1) + (mu1*(x2-x1)-mu1*(x-x1)) /2/pi/(x2-x1) * (theta2-theta1)
    wp_a = -(mu1*(x2-x1)-mu1*(x-x1)) /2/pi/(x2-x1) * log(r1/r2) + z/2/pi*(-mu1)/(x2-x1) * ((x2-x1)/z+theta2-theta1)

    # second panel point
    up_b = z/2/pi*(mu2)/(x2-x1)*log(r2/r1) + mu2*(x-x1)/2/pi/(x2-x1) * (theta2-theta1)
    wp_b = -(mu2*(x-x1))/2/pi/(x2-x1)*log(r1/r2) + z/2/pi*(mu2)/(x2-x1) * ((x2-x1)/z+theta2-theta1)

    vel_a = coordRot2D([up_a, wp_a], -theta, [0,0])
    vel_b = coordRot2D([up_b, wp_b], -theta, [0,0])

    vel = vel_a .+ vel_b
    
    return [vel; vel_a; vel_b]

end

pan_pts = readAF("airfoil.csv", true)
writeAF(pan_pts, "airfoil.dat")
pan_pts = repanel(pan_pts, 50, 1, true)
pan_pts, c_pts, thetas, norms, dists = procPanels(pan_pts)

linVort(1,1, [0,0], [1,0], [5,5])

# conditions
U = 1
chord = 1
alpha = 5
rho = 1.225

# Initialize solver matrix
num_pan = size(pan_pts,1) -1
A = zeros(num_pan+1, num_pan+1)
RHS = zeros(num_pan+1)
u_vec = U .* [cosd(alpha), sind(alpha)]

for i = 1:num_pan
    for j = 1:num_pan+1
        if i === j
            uw = linVort(1,1,pan_pts[j,1:2],pan_pts[j,3:4],c_pts[i])
        elseif j == num_pan+1
        else
        end
    end
end

# Kutta condition
A[end, 1] = 1
A[end, end] = 1
RHS[end] = 0