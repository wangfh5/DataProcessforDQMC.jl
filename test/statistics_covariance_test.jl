using Test
using DataProcessforDQMC: compute_cr_m2_covariance

@testset "compute_cr_m2_covariance" begin
    @testset "4-parameter version (original values)" begin
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

    @testset "3-parameter version (scaled values)" begin
        # Test 1: Basic calculation
        # scaled_m2=1.0, scaled_error=0.05, cr=0.5
        # σ_xy = 1.0 * (1-0.5) * (0.05/1.0)^2 = 0.00125
        result = compute_cr_m2_covariance(1.0, 0.05, 0.5)
        @test haskey(result, :covariance)
        @test !haskey(result, :scaled_m2)  # 3-param version doesn't return components
        @test isapprox(result.covariance, 0.00125; atol=1e-10)

        # Test 2: Zero error case
        result = compute_cr_m2_covariance(1.0, 0.0, 0.5)
        @test result.covariance == 0.0

        # Test 3: Edge case - cr=1
        result = compute_cr_m2_covariance(1.0, 0.05, 1.0)
        @test result.covariance == 0.0

        # Test 4: Edge case - cr=0
        result = compute_cr_m2_covariance(1.0, 0.05, 0.0)
        expected = 1.0 * 1.0 * (0.05/1.0)^2
        @test isapprox(result.covariance, expected; atol=1e-10)

        # Test 5: Input validation
        @test_throws AssertionError compute_cr_m2_covariance(-1.0, 0.05, 0.5)
        @test_throws AssertionError compute_cr_m2_covariance(1.0, -0.05, 0.5)
    end

    @testset "Equivalence between 3-param and 4-param versions" begin
        # The 3-param version with scaled values should give the same result
        # as the 4-param version with original values
        m2 = 2.0
        σ_m2 = 0.1
        cr = 0.5
        K = 0.5

        # 4-param version
        result_4p = compute_cr_m2_covariance(m2, σ_m2, cr, K)

        # 3-param version with pre-scaled values
        scaled_m2 = K * m2
        scaled_error = K * σ_m2
        result_3p = compute_cr_m2_covariance(scaled_m2, scaled_error, cr)

        @test isapprox(result_4p.covariance, result_3p.covariance; atol=1e-10)

        # Test with different parameters
        m2 = 5.0
        σ_m2 = 0.3
        cr = 0.7
        K = 0.2

        result_4p = compute_cr_m2_covariance(m2, σ_m2, cr, K)
        result_3p = compute_cr_m2_covariance(K * m2, K * σ_m2, cr)

        @test isapprox(result_4p.covariance, result_3p.covariance; atol=1e-10)
    end
end
