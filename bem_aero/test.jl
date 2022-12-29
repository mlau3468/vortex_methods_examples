include("bem_aero.jl")
using PAGE
using Plots

# a = vel_line_source_2d([0.6;1.0], [-1.0;3.0], [5.0;8.0])
# b = vel_line_source_2d_int([0.6;1.0], [-1.0;3.0], [5.0;8.0])
# a = pot_line_source_2d([0.0;0.0], [1.0;0.0], [2.0;2.0])
# b = pot_line_source_2d_int([0.0;0.0], [1.0;0.0], [2.0;2.0])
# a = pot_line_doublet_2d([1.0;2.0], [3.0;-1.0], [4.0;4.0])S
# b = pot_line_doublet_2d_int([1.0;2.0], [3.0;-1.0], [4.0;4.0])

# a = vel_line_doublet_2d([1.0;2.0], [3.0;-1.0], [4.0;4.0])
# b = vel_line_doublet_2d_int([1.0;2.0], [3.0;-1.0], [4.0;4.0])

# vela, velb = vel_line_vortex_linear_2d([1.0;2.0], [3.0;-1.0], [-4.0;4.0])
# vela2, velb2 = vel_line_vortex_linear_2d_int([1.0;2.0], [3.0;-1.0], [-4.0;4.0])

# phia, phib = pot_line_doublet_linear_2d([1.0;2.0], [3.0;-1.0], [-4.0;4.0])
# phia2, phib2 = pot_line_doublet_linear_2d_int([1.0;2.0], [3.0;-1.0], [-4.0;4.0])

phia, phib = pot_line_doublet_linear_2d([1.0;2.0], [3.0;-1.0], [-4.0;4.0])
phia2, phib2 = pot_line_doublet_linear_2d_int([1.0;2.0], [3.0;-1.0], [-4.0;4.0])
# phi = phia+phib
# a = pot_line_doublet_2d([1.0;2.0], [3.0;-1.0], [4.0;4.0])

phia, phib = pot_line_vortex_linear_2d([0;0], [1;1], [-2.0;2.0])
phia2, phib2 = pot_line_vortex_linear_2d_int([0;0], [1;1], [-2.0;2.0])

# quit()

phia, phib = pot_line_vortex_linear_2d([0;0], [1;0], [0.5;0])
phia2, phib2 = pot_line_vortex_linear_2d_self([0;0], [1;0], [0.5;0])

vela, velb = vel_line_vortex_linear_2d([2.0;1.0], [3.0;-1.0], [4.0;4.0])
vela2, velb2 = vel_line_vortex_linear_2d_ad([2.0;1.0], [3.0;-1.0], [4.0;4.0])

# phi1 = pot_line_source_2d([1.0;-2.0], [3.0;-1.0], [-4.0;4.0])
# phi2 = pot_line_source_2d_int([1.0;-2.0], [3.0;-1.0], [-4.0;4.0])

# phi1 = pot_line_vortex_2d([1.0;-2.0], [3.0;-1.0], [-4.0;4.0])
# phi2 = pot_line_vortex_2d_int([1.0;-2.0], [3.0;-1.0], [-4.0;4.0])

# phi = pot_line_vortex_2d([0;0], [2;2], [1;1])
# phi2 = pot_line_vortex_2d_self([0;0], [2;2], [1;1])

# phia, phib = pot_line_vortex_linear_2d([1;1], [2;2], [1;1])
# phia2, phib2 = pot_line_vortex_linear_2d_self([1;0], [2;2], [1;1])
# vela, velb = vel_line_vortex_linear_2d_self([0;0], [1;0], [0.5;0])
# vela2, velb2 = vel_line_vortex_linear_2d_self_ad([0;0], [1;0], [0.5;0])
