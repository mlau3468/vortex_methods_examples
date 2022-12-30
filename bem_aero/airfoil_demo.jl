include("bem_aero.jl")
using PAGE
using Plots

pan_vert = naca4(0.00, 0.0, 0.12, nchord=100, spacetype="cos", cosine_weight=1.0)
aoa = 2

result1 = airfoil_sourcedoublet_dirichlet(pan_vert, aoa)
result1 = airfoil_sourcedoublet_dirichlet(pan_vert, aoa, compute_full_potential=true)
result2 = airfoil_doublet_dirichlet(pan_vert, aoa)
result2 = airfoil_doublet_dirichlet(pan_vert, aoa, compute_full_potential=true)
result3 = airfoil_doublet_neumann(pan_vert, aoa)

result4 = airfoil_doublet_linear_dirichlet(pan_vert, aoa)
result5 = airfoil_vortex_linear_neumann(pan_vert, aoa)
result6 = airfoil_vortex_linear_neumann(pan_vert, aoa, num_integrate=true)

# result7 = airfoil_vortex_linear_neumann_pot(pan_vert, aoa)

cp_plot_const = plot(yflip=true, legend=:bottomright)
plot!(cp_plot_const, result1.pan_cpt[1,:], result1.pan_cp, label="source+doublet (dirichlet)")
plot!(cp_plot_const, result2.pan_cpt[1,:], result2.pan_cp, label="doublet (dirichlet)")
plot!(cp_plot_const, result3.pan_cpt[1,:], result3.pan_cp, label="doublet (neumann)")
xlabel!(cp_plot_const, "x")
ylabel!(cp_plot_const, "cp")
title!(cp_plot_const, "Constant Elements")

cp_plot_lin= plot(yflip=true, legend=:bottomright)
plot!(cp_plot_lin, result4.pan_cpt[1,:], result4.pan_cp, label="doublet (dirichlet)")
plot!(cp_plot_lin, result5.pan_cpt[1,:], result5.pan_cp, label="vortex (neumann)")
plot!(cp_plot_lin, result6.pan_cpt[1,:], result6.pan_cp, label="vortex (neumann, quadrature)")
# plot!(cp_plot_lin, result7.pan_cpt[1,:], result7.pan_cp, label="test")
xlabel!(cp_plot_lin, "x")
ylabel!(cp_plot_lin, "cp")
title!(cp_plot_lin, "Linear Elements")
# ylims!(cp_plot_lin, (-2, 2))

display(cp_plot_const)
display(cp_plot_lin)