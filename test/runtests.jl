using Test: @test, @testset

using ODE

@testset "ODE Tests" begin
    @testset "stepsize_adjust" begin
        (hnew, retry) = ODE.stepsize_adjust([1.0], [1e-8], 0.1, 1)
        @test isapprox(hnew, 0.1/(1e-8/(1e-8 + 1e-8))^0.5)
        @test retry == false

        (hnew, retry) = ODE.stepsize_adjust([2.234], [1e-3], 0.01, 3)
        @test isapprox(hnew, 0.5*0.01)
        @test retry == true

        (hnew, retry) = ODE.stepsize_adjust([2.234], [1e-12], 0.01, 3)
        @test isapprox(hnew, 2*0.01)
        @test retry == false
    end
end
