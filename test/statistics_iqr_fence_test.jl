using Test
using Statistics
using DataProcessforDQMC: iqr_fence_filter

@testset "iqr_fence_filter" begin
    # 1) Basic: remove one extreme outlier when n >= min_n
    values = [1.0, 2, 2, 2, 3, 3, 3, 1000]
    filtered, n_removed = iqr_fence_filter(values; k=10.0, min_n=8)
    @test n_removed == 1
    @test length(filtered) == 7
    @test maximum(filtered) < 100

    # 2) Not enough samples: do not apply fence
    values_small = [1.0, 2, 2, 3, 3, 3, 1000]
    filtered_small, n_removed_small = iqr_fence_filter(values_small; k=10.0, min_n=8)
    @test n_removed_small == 0
    @test maximum(filtered_small) == 1000.0

    # 3) IQR == 0: do not apply fence
    values_flat = fill(2.0, 8)
    filtered_flat, n_removed_flat = iqr_fence_filter(values_flat; k=10.0, min_n=8)
    @test n_removed_flat == 0
    @test filtered_flat == values_flat

    # 4) Ignore NaN/Inf/missing before filtering
    values_nf = Any[1.0, 2, 2, 2, 3, 3, 3, 4, NaN, Inf, missing, 1000]
    filtered_nf, n_removed_nf = iqr_fence_filter(values_nf; k=10.0, min_n=8)
    @test n_removed_nf == 1
    @test length(filtered_nf) == 8
    @test all(isfinite, filtered_nf)
    @test maximum(filtered_nf) <= 4.0
end
