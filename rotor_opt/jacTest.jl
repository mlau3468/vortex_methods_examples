using ReverseDiff: JacobianTape, JacobianConfig, jacobian, jacobian!, compile
using LinearAlgebra: mul!

#########
# setup #
#########
#=
# some objective functions to work with
f(a, b) = (a + b) * (a * b)'
# pre-record JacobianTapes for `f` and `g` using inputs of shape 10x10 with Float64 elements
const f_tape = JacobianTape(f, (rand(10, 10), rand(10, 10)))

# compile `f_tape` and `g_tape` into more optimized representations
const compiled_f_tape = compile(f_tape)

# some inputs and work buffers to play around with
a, b = rand(10, 10), rand(10, 10)
inputs = (a, b)
results = (similar(a, 100, 100), similar(b, 100, 100))
fcfg = JacobianConfig(inputs)

####################
# taking Jacobians #
####################

# with pre-recorded/compiled tapes (generated in the setup above) #
#-----------------------------------------------------------------#

# these should be the fastest methods, and non-allocating
jacobian!(results, compiled_f_tape, inputs)
=#

# some objective functions to work with
function f(a)
    return a^2
end

inputs = rand(10, 10)
results = similar(inputs, 100, 100)

# pre-record JacobianTapes for `f` and `g` using inputs of shape 10x10 with Float64 elements
const f_tape = JacobianTape(f, similar(inputs))

# compile `f_tape` and `g_tape` into more optimized representations
const compiled_f_tape = compile(f_tape)

# some inputs and work buffers to play around with

fcfg = JacobianConfig(inputs)

####################
# taking Jacobians #
####################

# with pre-recorded/compiled tapes (generated in the setup above) #
#-----------------------------------------------------------------#

# these should be the fastest methods, and non-allocating
jacobian!(results, compiled_f_tape, inputs)