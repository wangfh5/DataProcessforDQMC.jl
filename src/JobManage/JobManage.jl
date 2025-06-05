module JobManage

using Printf

#= 
Jobname contains information about what the Job is and all important paramaters, as precise as possible.
It is used as the directory name for the data of this job.
=#
include("JobNaming.jl")
export generate_jobname, parse_jobname, parse_jobname_legacy
export migrate_legacy_to_new, verify_migration

#=
A Job can have multiple simulations. These contain: 
- Different MPI processes in one run
- Another run with the same parameters
    - For update the data. 
    - For adding more measurements (typically Monte Carlo simulations)
=#
include("SimulationInfo.jl")
export SimulationInfo, save_simulation_info

end