include("singularity_elements.jl")

using Plots

a = pot_line_source_2d([0.6;1.0], [-1.0;3.0], [5.0;8.0])
b = pot_line_source_2d_int([0.6;1.0], [-1.0;3.0], [5.0;8.0])
