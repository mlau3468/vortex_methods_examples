using LinearAlgebra
using ForwardDiff
using FastGaussQuadrature

include("util.jl")
include("point_2d.jl")
include("line_2d.jl")
include("airfoil_panel_method/airfoil_panel_util.jl")
include("airfoil_panel_method/airfoil_source_doublet.jl")