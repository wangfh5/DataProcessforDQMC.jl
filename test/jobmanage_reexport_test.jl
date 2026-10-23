using DQMCJobManage: DQMCJobManage

@testset "DQMCJobManage re-export" begin
    @test DataProcessforDQMC.JobManage === DQMCJobManage
    for name in (:generate_jobname, :parse_jobname, :parse_jobname_legacy, :migrate_legacy_to_new,
                 :verify_migration, :SimulationInfo, :save_simulation_info)
        @test Base.isexported(DataProcessforDQMC, name)
        @test getfield(DataProcessforDQMC, name) === getfield(DQMCJobManage, name)
    end
    # transitional path used by old submit scripts and the sysimage precompile script
    @eval Main using DataProcessforDQMC.JobManage
    @test Main.generate_jobname("ft", [("prog", "submat", 0)]) == "ft_prog=submat"
end
