using DelimitedFiles
using CSV
include("aeroel.jl")
include("geo.jl")
include("vis.jl")
include("sim.jl")

U = 50
alpha = 5
uinf = [U*cos(deg2rad(alpha)); 0; U*sin(deg2rad(alpha))]

nspan = 13
nchord = 1
chord = 1

ys = [0;8]
chords = [chord;chord]
sweeps = [0;0]
twists =[0;0]
xref = 0.5

new_comp, new_pan = createWing("wing", ys, chords, sweeps, twists, nspan, nchord,xref,"cos")
panels = new_pan
panels2vtk(panels, "test.vtu")

c81 = readdlm("c81/naca0012.csv", ',', Float64)
