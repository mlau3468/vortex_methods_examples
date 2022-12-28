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
aoa = 0
result = airfoil_sourcedoublet_dirichlet(pan_vert, aoa)
result2 = airfoil_doublet_dirichlet(pan_vert, aoa)
result3 = airfoil_vortex_linear_neumann(pan_vert, aoa)

cp_plot = plot(yflip=true, legend=:bottomright)
plot!(cp_plot, result.pan_cpt[1,:], result.pan_cp, label="constant source+doublet")
plot!(cp_plot, result2.pan_cpt[1,:], result2.pan_cp, label="constant doublet")
plot!(cp_plot, result3.pan_cpt[1,:], result3.pan_cp, label="linear vortex")
xlabel!(cp_plot, "x")
ylabel!(cp_plot, "cp")