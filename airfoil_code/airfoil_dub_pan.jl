include("airfoil_util.jl")

function constDubPan(mu, p1, p2, p)
    #mu = doublet strength
    #p = point to calc velocity
    #p1 = start point of doublet panel
    #p2 = end point of doublet panel

    if p2 === nothing # wake
        theta = 0
        x1 = 0
        x = p[1] - p1[1]
        z = p[2] - p1[2]
        r_2 = (x-x1)^2+z^2
        u = mu/2/pi/r_2 * z
        w = -mu/2/pi / r_2 * x
    else
        theta = -atan(p2[2]-p1[2], p2[1]-p1[1]);
        p_new = coordRot2D(p, theta, p1)

        x1 = 0
        x2 = ptDist(p1, p2)

        x = p_new[1]
        z = p_new[2]
        if abs(z) < 0.00001 && x1<x && x<x2 # Influence on itself
            u = 0
            w = -mu/2/pi * (1/(x-x1) - 1/(x-x2))
        else
            u = mu/2/pi * (z/((x-x1)^2+z^2) - z/((x-x2)^2+z^2))
            w = -mu/2/pi * ((x-x1)/((x-x1)^2+z^2) - (x-x2)/((x-x2)^2+z^2))
        end
    end
    # transform back
    vel = [u, w]
    vel = coordRot2D(vel, -theta, [0,0])
    return vel

end

function airfoil_constDub(pan_pts, alpha)
pan_pts, c_pts, thetas, norms, tangents, dists = procPanels(pan_pts)
# append a panel for the wake
pan_pts = [pan_pts; [pan_pts[end,3], pan_pts[end,4], NaN, NaN]']

# conditions
U = 1
chord = 1
rho = 1.225

# Initialize solver matrix
num_pan = size(pan_pts,1) -1
A = zeros(num_pan+1, num_pan+1)
RHS = zeros(num_pan+1)
u_vec = U .* [cosd(alpha), sind(alpha)]

# Influence Matrix
for i in 1:num_pan
    RHS[i] =- u_vec'norms[i, :]
    for j in 1:num_pan+1
        if j == num_pan+1 # wake, last panel
            uw = constDubPan(1, pan_pts[j,1:2], nothing, c_pts[i,:])
        else
            uw = constDubPan(1, pan_pts[j,1:2], pan_pts[j,3:4], c_pts[i, :])
        end
        A[i,j] = uw'norms[i,:]
     end
end

# Kutta condition
A[end, 1] = 1
A[end, end] = 1
A[end, end-1] = -1
RHS[end] = 0

# remove one panel
remIdx = Int(round(size(pan_pts,1)/2))
num_pan = num_pan -1
keep_idx = [1:remIdx-1; remIdx+1:num_pan+2]

RHS = RHS[keep_idx]
A = A[keep_idx, keep_idx]
c_pts = c_pts[keep_idx[1:end-1], :]
norms = norms[keep_idx[1:end-1],:]
tangents = tangents[keep_idx[1:end-1],:]
thetas = thetas[keep_idx[1:end-1]]
pan_pts = pan_pts[keep_idx, :]
dists = dists[keep_idx[1:end-1]]

mu = A\RHS

pan_vels = zeros(num_pan, 2)
# velocities at each panel
for i in 1:num_pan
    pan_vels[i,:] = u_vec
    # perterbation velocity of panel on itself
    if i == 1  # First panel
        v_pan = 0.5*(mu[i+1]-mu[i])/ptDist(c_pts[i,:], c_pts[i+1,:])
    elseif i == num_pan  # last panel
        v_pan = 0.5*(mu[i]-mu[i-1])/ptDist(c_pts[i-1,:], c_pts[i,:])
    else
         v_pan = 0.5*(mu[i+1]-mu[i-1])/ptDist(c_pts[i-1,:], c_pts[i+1,:])
    end
    pan_vels[i, :] = pan_vels[i,:] +  coordRot2D([v_pan, 0], -thetas[i], [0,0])
    
    # influence of other panels
    for j in 1:num_pan+1
        if j == num_pan+1 # wake, last panel
            uw = constDubPan(mu[j], pan_pts[j,1:2], nothing, c_pts[i,:])
        else
            uw = constDubPan(mu[j], pan_pts[j,1:2], pan_pts[j,3:4], c_pts[i, :])
        end
        pan_vels[i, :] = pan_vels[i,:] + uw
     end
end

cp = 1 .- (pan_vels[:,1].^2 + pan_vels[:,2].^2) ./ U.^2
cpPlot = plot(c_pts[:,1], cp, yflip=true)

ps = cp.*0.5.*rho.*U.^2
fs = ps.*dists.*(-1).*norms[:,2]
cl = sum(fs)./0.5./chord./rho./U.^2
return [cl, cpPlot]
end


# Testing
pan_pts = readAF("airfoil.csv", true)
writeAF(pan_pts, "airfoil.dat")
#pan_pts = repanel(pan_pts, 50, 1, true)
aoa = 5
cl, cpPlot = airfoil_constDub(pan_pts, aoa)

println("CL: $cl")
display(cpPlot)