#import Pkg; Pkg.add("SafeTestsets")

using SafeTestsets
@safetestset "My f tests" begin include("my_f_test.jl") end
