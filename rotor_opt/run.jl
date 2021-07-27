using Base: Float64
using Plots: display
using DelimitedFiles
using Base: Float64
using ReverseDiff: GradientTape, GradientConfig, gradient, gradient!, compile, DiffResults, TrackedArray, TrackedReal
using ReverseDiff: JacobianTape, JacobianConfig, jacobian, jacobian!, compile
include("airfoilTools.jl")
include("ibl.jl")

#=
nacaNum = [2.0, 4.0, 1.0, 2.0]
pts = genNACA(nacaNum, 75)
res = airfoilCalc(pts, 30, 2)
display(res)
=#
function test(nacaNum)
    pts = genNACA(nacaNum, 75)
    res = airfoilCalc(pts, 500000, 2)
    res = [res[1] res[2]]
    return res
end

inputShape = [4, 1]
#outputShape = [(75)*2+1, 2]
outputShape = [2, 1]
gradShape = [prod(outputShape), prod(inputShape)]

inputs = [2.0, 4.0, 1.0, 2.0]
#display(test(inputs))
results = rand(gradShape...)

jacobian(test, inputs)

#=
# pre-record JacobianTapes for `f` and `g` using inputs of shape 10x10 with Float64 elements
const f_tape = JacobianTape(test, similar(inputs))

# compile `f_tape` and `g_tape` into more optimized representations
const compiled_f_tape = compile(f_tape)

# with pre-recorded/compiled tapes (generated in the setup above) #
#-----------------------------------------------------------------#

# these should be the fastest methods, and non-allocating
jac = jacobian!(results, compiled_f_tape, inputs)
display(jac)
=#