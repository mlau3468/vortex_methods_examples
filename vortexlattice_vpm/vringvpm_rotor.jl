using LinearAlgebra
using DelimitedFiles
include("aeroel.jl")
include("geo.jl")
include("vis.jl")
include("sim.jl")


nspan = 13
nchord = 4
ys = [1;2;4]
chords = [0.5;0.5;0.25]
sweeps = [0;0;10]
twists = 2 .*[12;8;4]
xref = 0.5

alpha = 0
U = 50
uinf = [U*cos(deg2rad(alpha)); 0; U*sin(deg2rad(alpha))]
uinf = [20;0;-5]
rho = 1.225
rpm = 500

om = rpm*2*pi/60
new_comp, new_pan = createWing("blade", ys, chords, sweeps, twists, nspan, nchord,xref)

new_comp.omega[:] = [0;0;om]


drev = 12
dt = 2*pi/om/drev
tsteps = drev*10
prefix = "test/_wing"



components = []
panels = []
te_idx = []

# add geometry
components = cat(components, new_comp, dims=1)
panels = cat(panels, new_pan, dims=1)
te_idx = cat(te_idx, new_comp.teidx, dims=1)


# ----------------------------------------------------------

# elements
particles = []
wakelines = []
wakerings = []

simulate(components, panels, te_idx, tsteps, dt, uinf, rho, particles, wakelines, wakerings, prefix)