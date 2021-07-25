using DelimitedFiles
using Plots
using Interpolations


function ptDist(pt1, pt2)
    return sqrt((pt2[1] - pt1[1])^2 + (pt2[2] - pt1[2])^2)
end

function coordRot2D(p, angle, origin)
    t = [cos(angle) -sin(angle); sin(angle) cos(angle)];
    return t * (p.-origin)
end

function readAF(filename, show=false)
    pts = readdlm(filename, ',', Float64)
    if show
        im = plot(pts[:,1], pts[:,2], aspect_ratio=1, line = (1, 2, :blue), marker = (:circle, 4, 0.6, :red));
        display(im)
    end
    num_pan = size(pts,1) -1
    pan_pts = zeros(num_pan, 4)
    pan_pts[1:num_pan, 1:2] = pts[1:num_pan,:]
    pan_pts[1:num_pan, 3:4] = pts[2:num_pan+1,:]
    return pan_pts
end

function writeAF(pan_pts, filename)
    open(filename, "w") do io
        writedlm(io, pan_pts[:, 1:2], ' ')
        writedlm(io, pan_pts[end,3:4]', ' ')
    end
end

function procPanels(pan_pts)
    pts = [pan_pts[:, 1:2];pan_pts[1, 1:2]']
    # collocation pts
    c_pts = (pts[1:end-1, :] + pts[2:end, :])./2;
    # normal vectors
    thetas = -atan.(pts[2:end, 2] - pts[1:end-1, 2], pts[2:end, 1] - pts[1:end-1, 1])
    norms = [sin.(thetas)'; cos.(thetas)']'
    tangents = [cos.(thetas)'; -sin.(thetas)']'
    #panel lengths
    dists = zeros(size(c_pts, 1))
    for i in 1:size(c_pts,1)
        dists[i] = ptDist(pts[i,:], pts[i+1,:])
    end
    return [pan_pts, c_pts, thetas, norms, tangents, dists]

end

function repanel(pan_pts, num_pan, wgt=0.5, show=false)
    pts = [pan_pts[:, 1:2];pan_pts[1, 1:2]']
    #pts = pan_pts[:, 1:2]
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
    if show
        im = plot(pts[:,1], pts[:,2], aspect_ratio=1, line = (1, 2, :blue), marker = (:circle, 4, 0.6, :red));
        display(im)
    end

    num_pan = size(pts,1) -1
    pan_pts = zeros(num_pan, 4)
    pan_pts[1:num_pan, 1:2] = pts[1:num_pan,:]
    pan_pts[1:num_pan, 3:4] = pts[2:num_pan+1,:]
    return pan_pts

end