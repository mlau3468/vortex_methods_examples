using ReverseDiff: JacobianTape, JacobianConfig, jacobian, jacobian!, compile

# some objective functions to work with
function f(a)
    return inv(a)
end

inputShape = [10 10]
outputShape = [10 10]
gradShape = [prod(outputShape), prod(outputShape)]

inputs = rand(inputShape...)
results = rand(gradShape...)

# pre-record JacobianTapes for `f` and `g` using inputs of shape 10x10 with Float64 elements
const f_tape = JacobianTape(f, similar(inputs))

# compile `f_tape` and `g_tape` into more optimized representations
const compiled_f_tape = compile(f_tape)

####################
# taking Jacobians #
####################

# with pre-recorded/compiled tapes (generated in the setup above) #
#-----------------------------------------------------------------#

# these should be the fastest methods, and non-allocating
jacobian!(results, compiled_f_tape, inputs)