mutable struct SimulationInfo
    # The path to the project directory, e.g., `~/Projects/MyProject`
    workdir::String
    # The path to the data directory, e.g., `~/Projects/MyProject/data`
    datadir::String
    # Set by JobManage.generate_jobname, contains the essential information about the job
    jobname::String
    # The path to the main directory of this job, e.g., `~/Projects/MyProject/data/jobname`
    maindir::String

    # process ID number (for MPI, label different processes in one run)
    MPIrank::Int
    # simulation ID number (for labling different identical runs; use it when you'd like to separate their results.)
    SimID::Int

    # Type of simulation: whether to resum or replace the old data (with the same SimID)
    # Typically, resum for Monte Carlo programs.
    # If there is no old data, it will be false. 
    resuming::Bool
end

function SimulationInfo(; workdir::String, jobname::String, datadir::String=joinpath(workdir, "data"), SimID::Int=0, MPIrank::Int=0, resuming::Bool=true)
    jobname = (SimID==0 ? jobname : join([jobname, "$(SimID)"], "_"))
    maindir = joinpath(datadir, jobname)

    if resuming # default setting
        resuming = isdir(maindir)
    else # if you insist on replacing the old data
        if isdir(maindir)
            @warn "The directory $maindir already exists and it will be replaced."
        end
    end

    return SimulationInfo(workdir, datadir, jobname, maindir, MPIrank, SimID, resuming)
end
function SimulationInfo(model::String, paras_info::Vector{Tuple{String,Real,Int}}; workdir::String, kwargs...)
    jobname = generate_jobname(model, paras_info)
    return SimulationInfo(workdir=workdir, jobname=jobname, kwargs...)
end

# print struct info as TOML format
function Base.show(io::IO, ::MIME"text/plain", sim_info::SimulationInfo)

    @printf io "[simulation_info]\n\n"
    @printf io "jobname  = \"%s\"\n" sim_info.jobname
    @printf io "SimID    = %d\n" sim_info.SimID
    @printf io "MPIrank  = %d\n" sim_info.MPIrank
    @printf io "julia_version     = \"%s\"" VERSION

    return nothing
end

function save_simulation_info(sim_info::SimulationInfo, additional_info = nothing)

    (; maindir, MPIrank, SimID) = sim_info
    fn = @sprintf "simulation_info_SimID%d_MPIrank%d.toml" SimID MPIrank
    open(joinpath(maindir, fn), "w") do fout
        show(fout, "text/plain", sim_info)
        if !isnothing(additional_info)
            @printf fout "\n\n"
            TOML.print(fout, Dict("additional_info" => additional_info), sorted = true)
        end
    end

    return nothing
end