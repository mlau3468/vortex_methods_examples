include("airfoil_util.jl")

function linVort(mu1, mu2, p1, p2, p)
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
    # first panel point
    up_a = z/2/pi*(-mu1)/(x2-x1)*log(r2/r1) + (mu1*(x2-x1)-mu1*(x-x1)) /2/pi/(x2-x1) * (theta2-theta1)
    wp_a = -(mu1*(x2-x1)-mu1*(x-x1)) /2/pi/(x2-x1) * log(r1/r2) + z/2/pi*(-mu1)/(x2-x1) * ((x2-x1)/z-(theta2-theta1))
    # katz plotkin pg 303. note error in wp term. a + is actually a -
    # second panel point
    up_b = z/2/pi*(mu2)/(x2-x1)*log(r2/r1) + mu2*(x-x1)/2/pi/(x2-x1) * (theta2-theta1)
    wp_b = -(mu2*(x-x1))/2/pi/(x2-x1)*log(r1/r2) + z/2/pi*(mu2)/(x2-x1) * ((x2-x1)/z-(theta2-theta1))

    vel_a = coordRot2D([up_a, wp_a], -theta, [0,0])
    vel_b = coordRot2D([up_b, wp_b], -theta, [0,0])

    vel = vel_a .+ vel_b
    
    return [vel'; vel_a'; vel_b']

end

pan_pts = readAF("airfoil.csv", true)
writeAF(pan_pts, "airfoil.dat")
#pan_pts = repanel(pan_pts, 50, 1.0, true)
pan_pts, c_pts, thetas, norms, tangents, dists = procPanels(pan_pts)
# conditions
U = 1
chord = 1
alpha = 5
rho = 1.225

# Initialize solver matrix
num_pan = size(pan_pts,1)
A = zeros(num_pan+1, num_pan+1) # normal matrix
B = zeros(num_pan, num_pan+1) # tangent matrix
RHS = zeros(num_pan+1)
u_vec = U .* [cosd(alpha), sind(alpha)]

for i = 1:num_pan
    RHS[i] = -u_vec'norms[i, :]
    for j = 1:num_pan + 1
        if j == 1
            uw = linVort(1,1,pan_pts[j,1:2],pan_pts[j,3:4],c_pts[i,:])
            uw_a = uw[2,:]
            uw_b = uw[3,:]
        elseif j == num_pan + 1
            uw = linVort(1,1,pan_pts[num_pan,1:2],pan_pts[num_pan,3:4],c_pts[i,:])
            uw_a = uw[2,:]
            uw_b = uw[3,:]
        else
            uw_b = linVort(1,1,pan_pts[j-1,1:2], pan_pts[j-1,3:4], c_pts[i,:])
            uw_b = uw_b[3,:]

            uw_a = linVort(1,1,pan_pts[j,1:2], pan_pts[j,3:4], c_pts[i,:])
            uw_a = uw_a[2,:]
        end
        uw = uw_a .+ uw_b
        A[i,j] = uw'norms[i,:]
        B[i,j] = uw'tangents[i,:]
    end
end

# fix b

for i = 1:num_pan
    for j = 1:num_pan
        mu1 = 1
        mu2 = 1
        p = c_pts[i,:]
        p1 = pan_pts[j,1:2]
        p2 = pan_pts[j,3:4]
        r1 = ptDist(p1, p)
        r2 = ptDist(p2, p)
        theta = -atan(p2[2]-p1[2], p2[1]-p1[1])
        p_new = coordRot2D(p, theta, p1)
    
        x1 = 0
        x2 = ptDist(p1, p2)
    
        
        x = p_new[1]
        z = p_new[2]

        if i == j
            up_a=-0.5*(x-x2)/(x2)
            up_b=0.5*(x)/(x2)
            wp_a=-0.15916
            wp_b=0.15916

            
            vel_a = coordRot2D([up_a, wp_a], -theta, [0,0])
            vel_b = coordRot2D([up_b, wp_b], -theta, [0,0])
        else
            theta1 = atan(z,x)
            theta2 = atan(z, x-x2)
            # first panel point
            up_a = z/2/pi*(-mu1)/(x2-x1)*log(r2/r1) + (mu1*(x2-x1)-mu1*(x-x1)) /2/pi/(x2-x1) * (theta2-theta1)
            wp_a = -(mu1*(x2-x1)-mu1*(x-x1)) /2/pi/(x2-x1) * log(r1/r2) + z/2/pi*(-mu1)/(x2-x1) * ((x2-x1)/z-(theta2-theta1))
            # katz plotkin pg 303. note error in wp term. a + is actually a -
            # second panel point
            up_b = z/2/pi*(mu2)/(x2-x1)*log(r2/r1) + mu2*(x-x1)/2/pi/(x2-x1) * (theta2-theta1)
            wp_b = -(mu2*(x-x1))/2/pi/(x2-x1)*log(r1/r2) + z/2/pi*(mu2)/(x2-x1) * ((x2-x1)/z-(theta2-theta1))
        
            vel_a = coordRot2D([up_a, wp_a], -theta, [0,0])
            vel_b = coordRot2D([up_b, wp_b], -theta, [0,0])
        end

        if j == 1
            #B[i,1] = vel_a'tangents[i,:]
            global holdb = vel_b'tangents[i,:]
        elseif j == num_pan
            #B[i,num_pan] = vel_a'tangents[i,:]+holdb
            #B[i,num_pan+1] = vel_b'tangents[i,:]
        else
            if i == j || i == j-1
            B[i,j] = vel_a'tangents[i,:]+holdb
            end
            global holdb=vel_b'tangents[i,:]
        end
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
println("CL: $cl")

# velocities
pan_vels = zeros(num_pan,2)
for i = 1:num_pan
    pan_vels[i,:] = u_vec
    #=
    temp = (mu[i] + mu[i+1])/4
    pan_vels[i,:] = pan_vels[i,:] + coordRot2D([temp, 0], thetas[i], [0,0])
    =#
    for j = 1:num_pan
        res = linVort(mu[j], mu[j+1], pan_pts[j,1:2], pan_pts[j,3:4], c_pts[i,:])
        pan_vels[i,:] = pan_vels[i,:] + res[1,:]
    end
end
cp = 1 .- (pan_vels[:,1].^2 + pan_vels[:,2].^2) ./ U.^2
im = plot(c_pts[:,1], cp, yflip=true)
display(im)

writedlm("B.csv", B, ',')

# fortran implementation
vs = zeros(num_pan)
cp = zeros(num_pan)
for i = 1:num_pan
    vs[i] = B[i,:]'mu + u_vec'tangents[i,:]
    cp[i] = 1-vs[i]^2/U^2
end
im = plot(c_pts[:,1], cp, yflip=true)
display(im)
