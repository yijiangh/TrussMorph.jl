using TrussMorph
using Test

# The testing setup
# Here we split different categories into tests
# Each category is a separate file
# By placing `@testset` on the file, it will run all of the tests there,
# and report back when the entire file is complete

# Optional: If you @testset begin ... end around the entire testing setup
# then all of your tests will run, and the failures will be counted
# and printed at the end. There is a balance between the completeness of
# the test results and the time it takes to run tests!

# Optional: tic() ... toc() around the tests as a quick way to keep track of
# large performance regressions.

@testset "TrussMorph.jl" begin
    include("parse_json_tests.jl")
    # include("stiffness_matrix_tests.jl")
    # include("viz_tests.jl")
end
