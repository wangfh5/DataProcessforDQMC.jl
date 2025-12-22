using Test
using DelimitedFiles
using DataProcessforDQMC: AFMCorrelationRatio

@testset "CorrelationRatio outlier filtering" begin
    mktempdir() do tmpdir
        data = zeros(Float64, 12, 4)  # kx, ky, real, imag

        # Q point bins
        for i in 1:6
            data[i, 1] = 0.0
            data[i, 2] = 0.0
        end
        data[1:6, 3] = [2.0, 2.1, 1.9, 2.0, 2.05, 100.0]

        # Q + delta q bins
        for i in 7:12
            data[i, 1] = 0.25
            data[i, 2] = 0.0
        end
        data[7:12, 3] = [1.0, 1.1, 0.9, 1.0, 1.05, 100.0]

        writedlm(joinpath(tmpdir, "afm_sf_k.bin"), data)

        # raw means (with outliers):
        # Q: (2.0+2.1+1.9+2.0+2.05+100)/6 = 110.05/6 ≈ 18.3416667
        # Q+δq: (1.0+1.1+0.9+1.0+1.05+100)/6 = 105.05/6 ≈ 17.5083333
        # R = 1 - 17.5083333/18.3416667 ≈ 0.0454
        res_raw = AFMCorrelationRatio((0.25, 0.0), (0.0, 0.0), "afm_sf_k.bin", tmpdir;
            startbin=1,
            outlier_mode=:dropmaxmin, outlier_param=0,
            verbose=false, auto_digits=false)
        @test isapprox(res_raw.S_AFM_Q, 110.05/6; atol=1e-6)
        @test isapprox(res_raw.S_AFM_Q_shifted, 105.05/6; atol=1e-6)
        @test isapprox(res_raw.correlation_ratio, 1 - (105.05/6)/(110.05/6); atol=1e-6)

        # iqrfence(k=1.5) removes only the 100.0 in each set:
        # Q mean = (2.0+2.1+1.9+2.0+2.05)/5 = 2.01
        # Q+δq mean = (1.0+1.1+0.9+1.0+1.05)/5 = 1.01
        # R ≈ 1 - 1.01/2.01 ≈ 0.4975
        res_iqr = AFMCorrelationRatio((0.25, 0.0), (0.0, 0.0), "afm_sf_k.bin", tmpdir;
            startbin=1,
            outlier_mode=:iqrfence, outlier_param=1.5,
            verbose=false, auto_digits=false)
        @test isapprox(res_iqr.S_AFM_Q, 2.01; atol=1e-6)
        @test isapprox(res_iqr.S_AFM_Q_shifted, 1.01; atol=1e-6)
        @test isapprox(res_iqr.correlation_ratio, 1 - 1.01/2.01; atol=1e-6)
    end
end
