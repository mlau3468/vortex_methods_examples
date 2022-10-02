import Base.@kwdef
using StaticArrays
using LinearAlgebra
using Interpolations

include(joinpath(@__DIR__, "src/geo.jl"))
include(joinpath(@__DIR__, "src/wake.jl"))
include(joinpath(@__DIR__, "src/pansys.jl"))
include(joinpath(@__DIR__, "src/aero.jl"))
include(joinpath(@__DIR__, "src/solver.jl"))
include(joinpath(@__DIR__, "src/vis.jl"))
