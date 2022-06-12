include("geo.jl")
include("aero.jl")
include("vis.jl")

xle = [0;0;0]
yle = [0;5;10]
zle = [0;0;0]
chords = [1;1;1]
pitch = [5;5;5]
nspan = [10;10]
tevec = [1;0;0]
telen = 30

xle, yle, zle, chord, pitch = refine_wing(xle, yle, zle, chords, pitch, nspan)

pansys = wing2liftline(xle, yle, zle, chord, pitch, tevec, telen)
vis_liftline(pansys.pan_con, pansys.pan_vert, "geometry")
alpha = deg2rad(0)
V = 1
vel_inf = [cos(alpha)*V;0;sin(alpha)*V]
CL, CDi = solve_liftline(pansys, vel_inf)
println("CL="*string(CL))
println("CDi="*string(CDi))

AR = 10
oswald = CL^2/(pi*AR*CDi)
println("e="*string(oswald))