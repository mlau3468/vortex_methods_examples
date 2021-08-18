using DelimitedFiles
using CSV
include("aeroel.jl")
include("geo.jl")
include("vis.jl")
include("sim.jl")

U = 50
alpha = 0
uinf = [U*cos(deg2rad(alpha)); 0; U*sin(deg2rad(alpha))]
rho = 1.225

nspan = 13
nchord = 1
chord = 1

ys = [0;6]
chords = [chord;chord]
sweeps = [0;0]
twists =[5;5]
xref = 0.5

new_comp, new_pan = createWing("wing", ys, chords, sweeps, twists, nspan, nchord,xref,"cos")
panels = new_pan
panels2vtk(panels, "test.vtu")

c81 = readdlm("c81/naca0012.csv", ',', Float64)

tsteps = 100
prefix = "test/_wing"
dt = chord/U/0.1


components = component[]
panels = vortRing[]
te_idx = []


# add geometry
components = cat(components, new_comp, dims=1)
panels = cat(panels, new_pan, dims=1)
te_idx = cat(te_idx, new_comp.teidx, dims=1)


# elements
particles = wakePart[]
wakelines = wakeLine[]
wakerings = wakeRing[]

#simulate(components, panels, te_idx, tsteps, dt, uinf, rho, particles, wakelines, wakerings, prefix)
simulate2(components, panels, te_idx, tsteps, dt, uinf, rho, particles, wakelines, wakerings, prefix)