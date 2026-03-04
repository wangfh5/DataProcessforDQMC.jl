using Test
using DelimitedFiles
using DataProcessforDQMC: StructureFactorAnalysis

@testset "StructureFactor periodic k-point matching" begin
    mktempdir() do tmpdir
        # Two k-points, each with >= 2 bins (rows).
        # k1 is on the BZ edge (0.5, 0.5). Shifting kx by +1/24 moves it beyond 0.5,
        # and the periodic-equivalent point should be (-0.5 + 1/24, 0.5) = (-0.45833333, 0.5).
        data = zeros(Float64, 4, 4)  # kx, ky, real, imag
        data[1:2, 1] .= 0.5
        data[1:2, 2] .= 0.5
        data[1:2, 3] .= 10.0
        data[1:2, 4] .= 0.0

        data[3:4, 1] .= -0.45833333
        data[3:4, 2] .= 0.5
        data[3:4, 3] .= 20.0
        data[3:4, 4] .= 0.0

        writedlm(joinpath(tmpdir, "spsm_k.bin"), data)

        res = StructureFactorAnalysis((0.5416666667, 0.5), "spsm_k.bin", tmpdir;
            startbin=1,
            outlier_mode=:dropmaxmin,
            outlier_param=0,
            auto_digits=false,
            verbose=false)

        @test isapprox(res.k_point[1], -0.45833333; atol=1e-6)
        @test isapprox(res.k_point[2], 0.5; atol=1e-6)
        @test isapprox(res.mean_real, 20.0; atol=1e-6)
    end
end

