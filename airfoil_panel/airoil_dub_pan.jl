using DelimitedFiles
using Plots
using Interpolations

function dub2D(mu, p1, p2, p)
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

function ptDist(pt1, pt2)
    return sqrt((pt2[1] - pt1[1])^2 + (pt2[2] - pt1[2])^2)
end

function coordRot2D(p, angle, origin)
    t = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    return t * (p.-origin)
end

function readAF(filename)
    pts = readdlm(filename, ',', Float64)
    num_pan = size(pts,1) -1
    pan_pts = zeros(num_pan + 1, 4)
    for i in 1:num_pan + 1
        if i == num_pan + 1
            pan_pts[i, :] = [pts[i, :]; NaN; NaN]
        else
            pan_pts[i, :] = [pts[i, :]; pts[i+1,:]]
        end
    end
    return pan_pts
end

function writeAF(pan_pts, filename)
    open(filename, "w") do io
        writedlm(io, pan_pts[:, 1:2], ' ')
    end
end

function procPanels(pan_pts)
    pts = pan_pts[:, 1:2]
    # collocation pts
    c_pts = (pts[1:end-1, :] + pts[2:end, :])./2;
    # normal vectors
    thetas = atan.(pts[2:end, 2] - pts[1:end-1, 2], pts[2:end, 1] - pts[1:end-1, 1])
    norms = [-sin.(thetas)'; cos.(thetas)']';
    #panel lengths
    dists = zeros(size(c_pts, 1))
    for i in 1:size(c_pts,1)
        dists[i] = ptDist(pts[i,:], pts[i+1,:])
    end
    return [pan_pts, c_pts, thetas, norms, dists]

end

function repanel(pan_pts, num_pan, wgt=0.5)
    pts = pan_pts[:, 1:2]
    dists = zeros(size(pts, 1)-1)
    for i in 1:size(pts,1)-1
        dists[i] = ptDist(pts[i,:], pts[i+1,:])
    end
    #find le point
    le_idx = findall(pts[:,1] .== minimum(pts[:,1]))
    if length(le_idx)>1
        le_idx = le_idx[round(length(le_idx)/2)]
    else
        le_idx = le_idx[1]
    end
    bot_pts = pts[1:le_idx,:]
    bot_dists = cumsum([0; dists[1:le_idx-1]])
    top_pts = pts[le_idx:end,:]
    top_dists = cumsum([0; dists[le_idx:end]])

    # distance weighting vector
    x = collect(range(0,stop=1,length=num_pan+1))
    x2 = (cos.(-pi.+pi.*x) .+ 1)./2
    x3 = x2.*wgt .+ x.*(1-wgt)
    
    #top points
    top_curve_x = LinearInterpolation(top_dists, top_pts[:,1])
    top_curve_y = LinearInterpolation(top_dists, top_pts[:,2])
    new_top_x = top_curve_x(top_dists[end] .* x3)
    new_top_y = top_curve_y(top_dists[end] .* x3)
    new_top_pts = hcat(new_top_x, new_top_y)

    #bottom points
    bot_curve_x = LinearInterpolation(bot_dists, bot_pts[:,1])
    bot_curve_y = LinearInterpolation(bot_dists, bot_pts[:,2])
    new_bot_x = bot_curve_x(bot_dists[end] .* x3)
    new_bot_y = bot_curve_y(bot_dists[end] .* x3)
    new_bot_pts = hcat(new_bot_x, new_bot_y)

    pts = [new_bot_pts; new_top_pts[2:end,:]]

    num_pan = size(pts,1) -1
    pan_pts = zeros(num_pan + 1, 4)
    for i in 1:num_pan + 1
        if i == num_pan + 1
            pan_pts[i, :] = [pts[i, :]; NaN; NaN]
        else
            pan_pts[i, :] = [pts[i, :]; pts[i+1,:]]
        end
    end
    return pan_pts

end

pan_pts = readAF("airfoil.csv")
writeAF(pan_pts, "airfoil.dat")
pan_pts = repanel(pan_pts, 50, 1)
pan_pts, c_pts, thetas, norms, dists = procPanels(pan_pts)

pts = pan_pts[:, 1:2]
im = plot(pts[:,1], pts[:,2], aspect_ratio=1, line = (1, 2, :blue), marker = (:circle, 4, 0.6, :red));
display(im)


# conditions
U = 1
chord = 1
alpha = 2
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
            uw = dub2D(1, pan_pts[j,1:2], nothing, c_pts[i,:])
        else
            uw = dub2D(1, pan_pts[j,1:2], pan_pts[j,3:4], c_pts[i, :])
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
remIdx = 30
num_pan = num_pan -1
keep_idx = [1:remIdx-1; remIdx+1:num_pan+2]

RHS = RHS[keep_idx]
A = A[keep_idx, keep_idx]
c_pts = c_pts[keep_idx[1:end-1], :]
norms = norms[keep_idx[1:end-1],:]
thetas = thetas[keep_idx[1:end-1]]
pan_pts = pan_pts[keep_idx, :]
dists = dists[keep_idx[1:end-1]]

mu = A\RHS

pan_vels = zeros(num_pan, 2)
# velocities at each panel
for i in 1:num_pan
    pan_vels[i,:] = u_vec
    if i == 1  # First panel
        v_pan = 0.5*(mu[i+1]-mu[i])/ptDist(c_pts[i,:], c_pts[i+1,:])
    elseif i == num_pan  # last panel
        v_pan = 0.5*(mu[i]-mu[i-1])/ptDist(c_pts[i-1,:], c_pts[i,:])
    else
         v_pan = 0.5*(mu[i+1]-mu[i-1])/ptDist(c_pts[i-1,:], c_pts[i+1,:])
    end
    pan_vels[i, :] = pan_vels[i,:] +  coordRot2D([v_pan, 0], thetas[i], [0,0])
    
    for j in 1:num_pan+1
        if j == num_pan+1 # wake, last panel
            uw = dub2D(mu[j], pan_pts[j,1:2], nothing, c_pts[i,:])
        else
            uw = dub2D(mu[j], pan_pts[j,1:2], pan_pts[j,3:4], c_pts[i, :])
        end
        pan_vels[i, :] = pan_vels[i,:] + uw
     end
end

cp = 1 .- (pan_vels[:,1].^2 + pan_vels[:,2].^2) ./ U.^2
ps = cp.*0.5.*rho.*U.^2
fs = ps.*dists.*(-1).*norms[:,2]
fs = sum(fs)./0.5./chord./rho./U.^2
println("CL: $fs")
plot(c_pts[:,1], cp, yflip=true)