using CauViz
using Test

@testset "CauViz.jl" begin
    @test my_f(2,1) == 5
    @test my_f(2,3) == 7
end

@testset "Derivative tests" begin
    @test derivative_of_my_f(2,1) == 2   
end
