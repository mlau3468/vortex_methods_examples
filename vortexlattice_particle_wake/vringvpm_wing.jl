using LinearAlgebra
using DelimitedFiles
include("aeroel.jl")
include("geo.jl")
include("vis.jl")
include("sim.jl")

nspan = 13
nchord = 4
chord = 1


ys = [0;8]
chords = [chord;chord]
sweeps = [0;0]
twists =[0;0]
xref = 0.5

new_comp, new_pan = createWing("wing", ys, chords, sweeps, twists, nspan, nchord,xref, "cos")
new_comp.vel[:] = [0;0;0]


alpha = 5
U = 50
uinf = [U*cos(deg2rad(alpha)); 0; U*sin(deg2rad(alpha))]
rho = 1.225

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

# test geometry motion
test_geo(components, panels, dt, prefix)

simulate(components, panels, te_idx, tsteps, dt, uinf, rho, particles, wakelines, wakerings, prefix)
#=
# lift coefficient
cl = cos(deg2rad(alpha))*total_force[3] - sin(deg2rad(alpha)) * total_force[1]
cl = cl/(1/2*rho*U^2*S)
println("Step: $t, CL=$cl")
=#