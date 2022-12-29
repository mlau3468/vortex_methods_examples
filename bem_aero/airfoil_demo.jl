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


cp_plot = plot(yflip=true, legend=:bottomright)
plot!(cp_plot, result1.pan_cpt[1,:], result1.pan_cp, label="constant source+doublet")
plot!(cp_plot, result2.pan_cpt[1,:], result2.pan_cp, label="constant doublet (dirichlet)")
plot!(cp_plot, result3.pan_cpt[1,:], result3.pan_cp, label="constant doublet (neumann)")

plot!(cp_plot, result4.pan_cpt[1,:], result4.pan_cp, label="linear doublet")
plot!(cp_plot, result5.pan_cpt[1,:], result5.pan_cp, label="linear vortex")
plot!(cp_plot, result6.pan_cpt[1,:], result6.pan_cp, label="linear vortex quadrature")


xlabel!(cp_plot, "x")
ylabel!(cp_plot, "cp")
# ylims!(cp_plot, (-2.0, 2.0))

display(cp_plot)

# plot(pan_vert[1,:], result5.vert_mu[1:end-1])