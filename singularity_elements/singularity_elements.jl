using LinearAlgebra
using ForwardDiff
using FastGaussQuadrature

include("elements_2d/util_2d.jl")
include("elements_2d/point_2d.jl")
include("elements_2d/line_2d.jl")

include("airfoil_panel_method/airfoil_panel_util.jl")
include("airfoil_panel_method/airfoil_source_doublet.jl")
include("airfoil_panel_method/airfoil_doublet.jl")
include("airfoil_panel_method/airfoil_vortex.jl")