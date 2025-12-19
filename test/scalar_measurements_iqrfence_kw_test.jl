using Test
using DataProcessforDQMC: EnergyAnalysis, RenyiNegativity
using DelimitedFiles

@testset "EnergyAnalysis outlier_mode/outlier_param" begin
    mktempdir() do tmpdir
        # Create a simple energy.bin with 10 rows, 9 columns
        # Normal values + one extreme outlier in the last row
        data = ones(10, 9)
        data[:, 5] .= [1.0, 1.1, 0.9, 1.0, 1.0, 1.1, 0.9, 1.0, 1.0, 1000.0]  # col 5 has outlier
        data[:, 7] .= 2.0
        data[:, 9] .= 3.0
        writedlm(joinpath(tmpdir, "energy.bin"), data)

        # Test dropmaxmin mode (trim 1 max and 1 min)
        result_trim = EnergyAnalysis("energy.bin", tmpdir;
            startbin=1, outlier_mode=:dropmaxmin, outlier_param=1, columns=[5], labels=["E"], verbose=false)
        # After trimming 1 max/min, outlier 1000 is removed
        @test result_trim["E"][1] < 100  # mean should be close to 1

        # Test iqrfence mode
        result_iqr = EnergyAnalysis("energy.bin", tmpdir;
            startbin=1, outlier_mode=:iqrfence, outlier_param=3.0, columns=[5], labels=["E"], verbose=false)
        # IQR fence should also remove the outlier
        @test result_iqr["E"][1] < 100

        # Test that dropmaxmin=0 (mode drop) without iqrfence keeps the outlier
        result_no_filter = EnergyAnalysis("energy.bin", tmpdir;
            startbin=1, outlier_mode=:dropmaxmin, outlier_param=0, columns=[5], labels=["E"], verbose=false)
        # Mean should be high due to the outlier
        @test result_no_filter["E"][1] > 50

        # Test endbin is inclusive
        result_endbin = EnergyAnalysis("energy.bin", tmpdir;
            startbin=1, endbin=5, outlier_mode=:dropmaxmin, outlier_param=0, columns=[5], labels=["E"], verbose=false)
        # Only first 5 rows, no outlier in those
        @test result_endbin["E"][1] < 2.0

        # Test negative dropmaxmin throws
        @test_throws ArgumentError EnergyAnalysis("energy.bin", tmpdir;
            startbin=1, outlier_mode=:dropmaxmin, outlier_param=-1, columns=[5], labels=["E"], verbose=false)
    end
end

@testset "RenyiNegativity outlier_mode/outlier_param" begin
    mktempdir() do tmpdir
        # Create expRenyiN3.bin with format: real imag real imag ... (L+1 pairs)
        # For L=2, we need 6 columns (3 pairs)
        # 10 bins, with variance in normal data and one outlier in the last row
        # Need variance so IQR > 0 for IQR fence to work
        data = ones(10, 6) * 0.5  # positive values for log
        data[:, 1] .= [0.4, 0.5, 0.6, 0.45, 0.55, 0.48, 0.52, 0.5, 0.5, 1000.0]  # real part col 1 has variance + outlier
        data[:, 3] .= 0.6
        data[:, 5] .= 0.7
        writedlm(joinpath(tmpdir, "expRenyiN3.bin"), data)

        # Test dropmaxmin mode
        exp_trim, renyi_trim = RenyiNegativity("expRenyiN3.bin", tmpdir;
            startbin=1, outlier_mode=:dropmaxmin, outlier_param=1, printLA=[])
        # First LA (column 1) should have mean < 100 after trimming
        @test exp_trim[1, 1] < 100

        # Test iqrfence mode - k=1.5 is a common mild outlier threshold
        exp_iqr, renyi_iqr = RenyiNegativity("expRenyiN3.bin", tmpdir;
            startbin=1, outlier_mode=:iqrfence, outlier_param=1.5, printLA=[])
        @test exp_iqr[1, 1] < 100

        # Test negative dropmaxmin throws
        @test_throws ArgumentError RenyiNegativity("expRenyiN3.bin", tmpdir;
            startbin=1, outlier_mode=:dropmaxmin, outlier_param=-1)
    end
end
