using Test
using DelimitedFiles
using DataProcessforDQMC: CDWStructureFactor, create_cdw_from_nn

@testset "CDWStructureFactor regenerates a stale cdw_k.bin" begin
    mktempdir() do tmpdir
        # one k-point, columns: kx, ky, then (re, im) of AA, BA, AB, BB
        row(aa) = [0.0 0.0 aa 0.0 -1.0 0.0 -1.0 0.0 aa 0.0]
        nn_path = joinpath(tmpdir, "nn_k.bin")

        # cdw_k.bin derived when the run had two bins ...
        writedlm(nn_path, vcat(row(1.0), row(2.0)))
        create_cdw_from_nn("nn_k.bin", "cdw_k.bin", tmpdir, tmpdir; verbose=false)
        # ... then the run continued and appended two more bins
        writedlm(nn_path, vcat(row(1.0), row(2.0), row(3.0), row(4.0)))

        res = CDWStructureFactor((0.0, 0.0), "cdw_sf_k.bin", tmpdir;
            source_file="cdw_k.bin", force_rebuild=true, startbin=1,
            outlier_mode=:dropmaxmin, outlier_param=0, auto_digits=false, verbose=false)

        @test countlines(joinpath(tmpdir, "cdw_k.bin")) == 4
        @test countlines(joinpath(tmpdir, "cdw_sf_k.bin")) == 4
        # S = AA + BB - BA - AB = 2aa + 2, averaged over all four bins
        @test isapprox(res.mean_real, 7.0; atol=1e-6)
    end
end
