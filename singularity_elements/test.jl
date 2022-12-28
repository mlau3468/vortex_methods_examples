include("singularity_elements.jl")

using Plots

# a = vel_line_source_2d([0.6;1.0], [-1.0;3.0], [5.0;8.0])
# b = vel_line_source_2d_int([0.6;1.0], [-1.0;3.0], [5.0;8.0])
# a = pot_line_source_2d([0.0;0.0], [1.0;0.0], [2.0;2.0])
# b = pot_line_source_2d_int([0.0;0.0], [1.0;0.0], [2.0;2.0])
a = vel_line_vortex_2d([0.0;0.0], [1.0;0.0], [2.0;2.0])
b = vel_line_vortex_2d_int([0.0;0.0], [1.0;0.0], [2.0;2.0])
