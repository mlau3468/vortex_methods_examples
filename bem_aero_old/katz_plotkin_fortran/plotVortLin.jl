using Plots

pts = readdlm("CPLV.DAT", ',', Float64)
plot(pts[:,1], pts[:,2], yflip=true)