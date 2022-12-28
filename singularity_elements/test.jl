include("singularity_elements.jl")

using Plots

a = pot_line_source_2d([0;1], [1;0], [2;2])
b = pot_line_source_2d_int([0;1], [1;0], [2;2])
