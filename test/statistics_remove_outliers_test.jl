using Test
using DataProcessforDQMC: remove_outliers

@testset "remove_outliers" begin
    @testset "dropmaxmin mode" begin
        vals = [1.0, 2, 3, 4, 100.0]
        res = remove_outliers(vals, :dropmaxmin, 1; min_n=5)
        @test res.values == [2.0, 3.0, 4.0]
        @test res.removed_min == 1
        @test res.removed_max == 1

        # min_n gating: too short, no filtering
        short_vals = [1.0, 100.0, 2.0]
        res_short = remove_outliers(short_vals, :dropmaxmin, 1; min_n=5)
        @test res_short.values == short_vals
        @test res_short.removed_min == 0
        @test res_short.removed_max == 0
    end

    @testset "iqrfence mode" begin
        vals = [1.0, 2, 2, 2, 3, 3, 3, 1000.0]
        res = remove_outliers(vals, :iqrfence, 10.0; min_n=5)
        @test maximum(res.values) < 100.0
        @test length(res.values) == 7
        @test res.removed_min == 0
        @test res.removed_max == 1

        # min_n gating
        short_vals = [1.0, 2, 100.0]
        res_short = remove_outliers(short_vals, :iqrfence, 10.0; min_n=5)
        @test res_short.values == short_vals
        @test res_short.removed_min == 0
        @test res_short.removed_max == 0
    end

    @test_throws ArgumentError remove_outliers([1,2,3], :iqrfence, -1.0)
    @test_throws ArgumentError remove_outliers([1,2,3], :dropmaxmin, -1)
    @test_throws ArgumentError remove_outliers([1,2,3], :unknown, 1)
end
