include("singularity_elements.jl")
using PAGE
using Plots

# a = vel_line_source_2d([0.6;1.0], [-1.0;3.0], [5.0;8.0])
# b = vel_line_source_2d_int([0.6;1.0], [-1.0;3.0], [5.0;8.0])
# a = pot_line_source_2d([0.0;0.0], [1.0;0.0], [2.0;2.0])
# b = pot_line_source_2d_int([0.0;0.0], [1.0;0.0], [2.0;2.0])
a = pot_line_doublet_2d([1.0;2.0], [3.0;-1.0], [4.0;4.0])
b = pot_line_doublet_2d_int([1.0;2.0], [3.0;-1.0], [4.0;4.0])

a = vel_line_doublet_2d([1.0;2.0], [3.0;-1.0], [4.0;4.0])
b = vel_line_doublet_2d_int([1.0;2.0], [3.0;-1.0], [4.0;4.0])

pan_vert = naca4(0.00, 0.0, 0.04, nchord=100, spacetype="cos", cosine_weight=1.0)
aoa = 1
result = airfoil_sourcedoublet_dirichlet(pan_vert, aoa)
result2 = airfoil_doublet_dirichlet(pan_vert, aoa)
# result3 = airfoil_vortex_neumann(pan_vert, aoa)

cp_plot = plot(yflip=true, legend=:bottomright)
plot!(cp_plot, result.pan_cpt[1,:], result.pan_cp, label="sourc+doub (dir)")
plot!(cp_plot, result2.pan_cpt[1,:], result2.pan_cp, label="doub (dir)")
# plot!(cp_plot, result3.pan_cpt[1,:], result3.pan_cp, label="doub (neu)")
xlabel!(cp_plot, "x")
ylabel!(cp_plot, "cp")

# function coordRot2D(p, angle, origin)
#     t = [cos(angle) -sin(angle); sin(angle) cos(angle)];
#     return t * (p.-origin)
# end

# function ptDist(pt1, pt2)
#     return sqrt((pt2[1] - pt1[1])^2 + (pt2[2] - pt1[2])^2)
# end

# function constDubPan(mu, p1, p2, p)
#     #mu = doublet strength
#     #p = point to calc velocity
#     #p1 = start point of doublet panel
#     #p2 = end point of doublet panel

#     if p2 === nothing # wake
#         theta = 0
#         x1 = 0
#         x = p[1] - p1[1]
#         z = p[2] - p1[2]
#         r_2 = (x-x1)^2+z^2
#         u = mu/2/pi/r_2 * z
#         w = -mu/2/pi / r_2 * x
#     else
#         theta = -atan(p2[2]-p1[2], p2[1]-p1[1]);
#         p_new = coordRot2D(p, theta, p1)

#         x1 = 0
#         x2 = ptDist(p1, p2)

#         x = p_new[1]
#         z = p_new[2]
#         if abs(z) < 0.00001 && x1<x && x<x2 # Influence on itself
#             u = 0
#             w = -mu/2/pi * (1/(x-x1) - 1/(x-x2))
#         else
#             u = mu/2/pi * (z/((x-x1)^2+z^2) - z/((x-x2)^2+z^2))
#             w = -mu/2/pi * ((x-x1)/((x-x1)^2+z^2) - (x-x2)/((x-x2)^2+z^2))
#         end
#     end
#     # transform back
#     vel = [u, w]
#     vel = coordRot2D(vel, -theta, [0,0])
#     return vel

# end

# a = constDubPan(1, [0;0], [1;0], [3;2])
# b = vel_line_doublet_2d([0;0], [1;0], [3;2])