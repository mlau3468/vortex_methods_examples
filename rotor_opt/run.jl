using Base: Float64
using Plots: display
using DelimitedFiles
using Base: Float64
using ReverseDiff: GradientTape, GradientConfig, gradient, gradient!, compile, DiffResults, TrackedArray, TrackedReal
include("airfoilTools.jl")
include("ibl.jl")


function test(nacaNum)
    pts = genNACA(nacaNum, 75)
    return sum(pts)
end


inputs = [2.0,4.0,1.0,2.0]

# pre-record a GradientTape for `f` using inputs of shape 100x100 with Float64 elements
const f_tape = GradientTape(test, (similar(inputs)))

# compile `f_tape` into a more optimized representation
const compiled_f_tape = compile(f_tape)

# some inputs and work buffers to play around with

results = similar(inputs)
all_results = map(DiffResults.GradientResult, [results])
cfg = GradientConfig(inputs)

####################
# taking gradients #
####################

# with pre-recorded/compiled tapes (generated in the setup above) #
#-----------------------------------------------------------------#

# this should be the fastest method, and non-allocating
grad = gradient!(results, compiled_f_tape, inputs)
display(grad)
