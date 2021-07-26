using Base: Float64
using Plots: display
using DelimitedFiles
using ForwardDiff: jacobian, jacobian!, JacobianConfig, Chunk
using BenchmarkTools
include("airfoilTools.jl")
include("ibl.jl")

function test!(output,nacaNum)
    pts = genNACA(nacaNum, 75)
    res = airfoilCalc(pts, 30, 2)
    output[1] = res[1]
    output[2] = res[2]
    return output
end

inputShape = [4, 1]
outputShape = [2, 1]
gradShape = [prod(outputShape), prod(inputShape)]

inputs = [2.0, 4.0, 1.0, 2.0]
outputs = Array{Float64}(undef, outputShape...)
gradresults = rand(gradShape...)


const jconfig = JacobianConfig(test!, outputs, inputs, Chunk{4}())

jacobian!(gradresults,test!, outputs,inputs, jconfig)

display(outputs)
display(gradresults)