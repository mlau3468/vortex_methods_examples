using Base: Float64
using ReverseDiff: GradientTape, GradientConfig, gradient, gradient!, compile, DiffResults, TrackedArray, TrackedReal


#=
# some objective function to work with
f(a, b) = sum(a' * b + a * b')

# pre-record a GradientTape for `f` using inputs of shape 100x100 with Float64 elements
const f_tape = GradientTape(f, (rand(100, 100), rand(100, 100)))

# compile `f_tape` into a more optimized representation
const compiled_f_tape = compile(f_tape)

# some inputs and work buffers to play around with
a, b = rand(100, 100), rand(100, 100)
inputs = (a, b)
results = (similar(a), similar(b))
all_results = map(DiffResults.GradientResult, results)
cfg = GradientConfig(inputs)

####################
# taking gradients #
####################

# with pre-recorded/compiled tapes (generated in the setup above) #
#-----------------------------------------------------------------#

# this should be the fastest method, and non-allocating
grad = gradient!(results, compiled_f_tape, inputs)
display(grad)
=#

# some objective function to work with
function f(a)
    #test = Array{TrackedReal}(undef,size(a))
    test = similar(a, size(a))
    for i = 1:size(a,1)
        for j =1:size(a,2)
            test[i,j] = a[i,j]
        end
    end
    test[1,1] = 1
    return sum(test.^2)
end

inputs = rand(10, 10)
results = similar(inputs)

# pre-record a GradientTape for `f` using inputs of shape 100x100 with Float64 elements
const f_tape = GradientTape(f, (similar(inputs)))

# compile `f_tape` into a more optimized representation
const compiled_f_tape = compile(f_tape)

# some inputs and work buffers to play around with


all_results = map(DiffResults.GradientResult, [results])
cfg = GradientConfig(inputs)

# Take gradients with pre-recorded/compiled tapes (generated in the setup above) #
#-----------------------------------------------------------------#

# this should be the fastest method, and non-allocating
grad = gradient!(results, compiled_f_tape, inputs)
display(inputs)
display(grad)
