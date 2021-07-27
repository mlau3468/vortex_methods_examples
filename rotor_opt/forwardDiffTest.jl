using Base: Float64
using Plots: display
using DelimitedFiles
using ForwardDiff: jacobian, jacobian!, JacobianConfig, Chunk
using BenchmarkTools
include("airfoilTools.jl")
include("ibl.jl")

function test!(output,nacaNum)
    pts = genNACA(nacaNum, 75)
    alphas = 1:2
    for a = 1:length(alphas)
        res = airfoilCalc(pts, 30, a)
        output[a,1] = res[1]
        output[a,2] = res[2]
    end
    return output
end

inputShape = [4, 1]
outputShape = [2, 2]
gradShape = [prod(outputShape), prod(inputShape)]

inputs = [2.0, 4.0, 1.0, 2.0]
outputs = similar(inputs, outputShape...)
gradresults = similar(inputs, gradShape...)


const jconfig = JacobianConfig(test!, outputs, inputs, Chunk{4}())

jacobian!(gradresults,test!, outputs,inputs, jconfig)
display(outputs)
display(gradresults)