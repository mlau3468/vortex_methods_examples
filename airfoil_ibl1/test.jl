using Plots: display
using DelimitedFiles
include("airfoilTools.jl")
include("ibl.jl")

pts = genNACA([2,4,1,2], 75)
p = plot(pts[:,1], pts[:,2], aspect_ratio=:equal)
display(p)

cl, cd, err = airfoilCalc(pts, 30, 10)
println(cl)
println(cd)