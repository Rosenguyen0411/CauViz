using CauViz
using Test

@testset "CauViz.jl" begin
    @test my_f(2,1) == 5
    @test my_f(2,3) == 7
    @test my_f(2,3) == 19
end
