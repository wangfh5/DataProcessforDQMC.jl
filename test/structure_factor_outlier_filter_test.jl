using Test
using DelimitedFiles
using DataProcessforDQMC: StructureFactorAnalysis

@testset "StructureFactor outlier filtering" begin
    mktempdir() do tmpdir
        data = zeros(Float64, 6, 4)  # kx, ky, real, imag
        data[:, 1] .= 0.0
        data[:, 2] .= 0.0
        data[:, 3] .= [1.0, 1.1, 0.9, 1.0, 1.0, 50.0]  # one extreme outlier
        data[:, 4] .= 0.0
        writedlm(joinpath(tmpdir, "spsm_k.bin"), data)

        # baseline: no filtering, mean dominated by the outlier (expect 55/6)
        res_raw = StructureFactorAnalysis((0.0, 0.0), "spsm_k.bin", tmpdir;
            startbin=1, outlier_mode=:dropmaxmin, outlier_param=0, auto_digits=false, verbose=false)
        @test isapprox(res_raw.mean_real, 55/6; atol=1e-6)

        # dropmaxmin removes one min (0.9) and one max (50.0) -> mean = (1+1+1+1.1)/4 = 1.025
        res_trim = StructureFactorAnalysis((0.0, 0.0), "spsm_k.bin", tmpdir;
            startbin=1, outlier_mode=:dropmaxmin, outlier_param=1, auto_digits=false, verbose=false)
        @test isapprox(res_trim.mean_real, 1.025; atol=1e-6)

        # iqrfence(k=1.5) removes only the 50.0 outlier -> mean = 5.0/5 = 1.0
        res_iqr = StructureFactorAnalysis((0.0, 0.0), "spsm_k.bin", tmpdir;
            startbin=1, outlier_mode=:iqrfence, outlier_param=1.5, auto_digits=false, verbose=false)
        @test isapprox(res_iqr.mean_real, 1.0; atol=1e-6)
    end
end
