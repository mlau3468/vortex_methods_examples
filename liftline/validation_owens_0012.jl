using Interpolations
using LinearAlgebra
using DelimitedFiles
using XLSX
using DataFrames
using Plots
include("geo.jl")
include("aero.jl")
include("vis.jl")

# Create Geometry
span = 232.02
xle = [0;0;0].*0.0254
yle = [-span/2;0;span/2].*0.0254
zle = [0;0;0].*0.0254
chords = [39.37;39.37;39.37].*0.0254
area = 9134.6*0.00064516
AR = 5.8933
pitch = [0;0;0]
nspan = 50
telen = 30*sum(chords)./length(chords)
xle, yle, zle, chord, pitch = refine_wing(xle, yle, zle, chords, pitch, nspan, "cos")

# Freestream
alphas = collect(LinRange(-4,20,30))
V = 10
cls = zeros(length(alphas))
cds = zeros(length(alphas))

for i = 1:length(alphas)
    # Freestream definition
    alpha = deg2rad(alphas[i])
    vel_inf = [cos(alpha)*V;0;sin(alpha)*V]
    v_mag = norm(vel_inf)
    rho = 1.225
    tevec = [cos(alpha);0;sin(alpha);]
    pansys = wing2liftline(xle, yle, zle, chord, pitch, tevec, telen)
    vis_liftline(pansys.pan_con, pansys.pan_vert, "geometry")

    # Solve lifting line, Wickenheiser
    #F_inv, F_visc = solve_liftline_weissinger(pansys, vel_inf, rho, "naca0012_mod.csv")
    F_inv, F_visc = solve_liftline_vandam(pansys, vel_inf, rho, "naca0012_mod.csv")
    lift = F_visc[3]*cos(alpha) - F_visc[1]*sin(alpha)
    drag = F_visc[3]*sin(alpha) + F_visc[1]*cos(alpha)
    CL = lift/(1/2*rho*v_mag^2*area)
    CD = drag/(1/2*rho*v_mag^2*area)
    cls[i] = CL
    cds[i] = CD
end

cds .= cds .+ 0.003

# windtunnel data
df = DataFrame(XLSX.readtable("validation_owens_0012.xlsx", "Sheet1")...)

p1 = plot(alphas, cls, label="model", legend=:bottomright)
scatter!(df[!,"alpha"], df[!,"cl"], label="data")
xlabel!("Angle of Attack (deg)")
ylabel!("CL")
ylims!((-0.4,1.2))
display(p1)

println(maximum(cls))
println(maximum(df[!,"cl"]))

p2 = plot(cds, cls, label="model", legend=:bottomright)
scatter!(df[!,"cd_polar"], df[!,"cl_polar"], label="data")
xlabel!("CD")
ylabel!("CL")
display(p2)