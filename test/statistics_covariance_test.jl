using Test
using DataProcessforDQMC: compute_cr_m2_covariance

@testset "compute_cr_m2_covariance" begin
    # Test 1: Basic calculation + return_components flag
    # m2=2.0, σ_m2=0.1, cr=0.5, K=0.5
    # σ_xy = 1.0 * (1-0.5) * (0.1/2.0)^2 = 0.00125
    result = compute_cr_m2_covariance(2.0, 0.1, 0.5, 0.5; return_components=true)
    @test haskey(result, :covariance)
    @test haskey(result, :scaled_m2)
    @test haskey(result, :relative_error_m2)
    @test isapprox(result.covariance, 0.00125; atol=1e-10)
    @test isapprox(result.scaled_m2, 1.0; atol=1e-10)
    @test isapprox(result.relative_error_m2, 0.05; atol=1e-10)

    # Test 2: Zero error case
    result = compute_cr_m2_covariance(2.0, 0.0, 0.5, 0.5)
    @test result.covariance == 0.0

    # Test 3: Edge case - cr=1 (perfect correlation, B=0)
    # When cr=1, x=1, so (1-x)=0, thus covariance=0
    result = compute_cr_m2_covariance(2.0, 0.1, 1.0, 0.5)
    @test result.covariance == 0.0

    # Test 4: Edge case - cr=0 (no correlation, B=A)
    # When cr=0, x=0, so (1-x)=1
    # σ_xy = y * 1.0 * (σ_A/A)²
    result = compute_cr_m2_covariance(2.0, 0.1, 0.0, 0.5)
    expected = 0.5 * 2.0 * 1.0 * (0.1/2.0)^2
    @test isapprox(result.covariance, expected; atol=1e-10)

    # Test 5: Input validation
    @test_throws AssertionError compute_cr_m2_covariance(-1.0, 0.1, 0.5, 0.5)
    @test_throws AssertionError compute_cr_m2_covariance(2.0, -0.1, 0.5, 0.5)
    @test_throws AssertionError compute_cr_m2_covariance(2.0, 0.1, 0.5, 0.0)
end
