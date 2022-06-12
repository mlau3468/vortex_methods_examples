using Interpolations
using LinearAlgebra
include("geo.jl")
include("aero.jl")
include("vis.jl")

# Freestream definition
alpha = deg2rad(5)
V = 1
vel_inf = [cos(alpha)*V;0;sin(alpha)*V]
v_mag = norm(vel_inf)
rho = 1.225

# Create Geometry
xle = [0;0;0]
yle = [0;5;10]
zle = [0;0;0]
chords = [1;1;1]
area = 10
AR = 10
pitch = [0;0;0]
nspan = [10;10]
tevec = [cos(alpha);0;sin(alpha);]
telen = 30*sum(chords)./length(chords)

xle, yle, zle, chord, pitch = refine_wing(xle, yle, zle, chords, pitch, nspan)
pansys = wing2liftline(xle, yle, zle, chord, pitch, tevec, telen)
vis_liftline(pansys.pan_con, pansys.pan_vert, "geometry")

# Solve lifting line, Katz Plotkin
lift, drag = solve_liftline(pansys, vel_inf, rho)
CL = lift/(1/2*rho*v_mag^2*area)
CDi = drag/(1/2*rho*v_mag^2*area)
println("Katz-Plotkin Method")
println("CL="*string(CL))
println("CDi="*string(CDi))
println("")

# Solve lifting line, modified Katz Plotkin
F = solve_liftline2(pansys, vel_inf, rho)
lift = F[3]*cos(alpha) - F[1]*sin(alpha)
drag = F[3]*sin(alpha) + F[1]*cos(alpha)
CL = lift/(1/2*rho*v_mag^2*area)
CDi = drag/(1/2*rho*v_mag^2*area)
println("Modified Katz-Plotkin Method")
println("CL="*string(CL))
println("CDi="*string(CDi))
println("")

oswald = CL^2/(pi*AR*CDi)
println("e="*string(oswald))