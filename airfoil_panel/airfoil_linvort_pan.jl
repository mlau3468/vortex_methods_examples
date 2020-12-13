include("airfoil_util.jl")

function linVort(mu1, mu2, p1, p2, p)
    
end

pan_pts = readAF("airfoil.csv", true)
writeAF(pan_pts, "airfoil.dat")
pan_pts = repanel(pan_pts, 50, 1, true)
pan_pts, c_pts, thetas, norms, dists = procPanels(pan_pts)
