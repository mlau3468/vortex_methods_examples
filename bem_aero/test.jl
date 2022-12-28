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

vela, velb = vel_line_vortex_linear_2d([1.0;2.0], [3.0;-1.0], [4.0;4.0])
vel = vela.*0 .+ velb.*1
vel2 = vel_line_vortex_linear_2d_int([1.0;2.0], [3.0;-1.0], [4.0;4.0])