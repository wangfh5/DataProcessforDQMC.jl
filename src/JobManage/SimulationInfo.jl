"""
    mutable struct SimulationInfo

A mutable struct to store simulation information including paths, job details, and execution context.

# Fields
- `workdir::String`: The path to the project directory, e.g., `~/Projects/MyProject`
- `datadir::String`: The path to the data directory, e.g., `~/Projects/MyProject/data`
- `jobname::String`: Set by JobManage.generate_jobname, contains the essential information about the job
- `maindir::String`: The path to the main directory of this job, e.g., `~/Projects/MyProject/data/jobname`
- `MPIrank::Int`: Process ID number (for MPI, label different processes in one run)
- `SimID::Int`: Simulation ID number (for labeling different identical runs; use it when you'd like to separate their results)
- `resuming::Bool`: Type of simulation: whether to resume or replace the old data (with the same SimID)
"""
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

"""
    SimulationInfo(; workdir::AbstractString, jobname::AbstractString, datadir::AbstractString=joinpath(workdir, "data"), SimID::Int=0, MPIrank::Int=0, resuming::Bool=true)

Create a SimulationInfo object with keyword arguments.

# Arguments
- `workdir::AbstractString`: The path to the project directory
- `jobname::AbstractString`: The jobname containing essential job information
- `datadir::AbstractString`: The path to the data directory (default: `joinpath(workdir, "data")`)
- `SimID::Int`: Simulation ID number (default: 0)
- `MPIrank::Int`: Process ID number for MPI (default: 0)
- `resuming::Bool`: Whether to resume existing simulation (default: true)

# Returns
- `SimulationInfo`: The constructed simulation info object

# Examples
```julia
sim_info = SimulationInfo(
    workdir="/path/to/project",
    jobname="honeycomb_L=8_U=4.0_enable_afm=T"
)
```
"""
function SimulationInfo(; workdir::AbstractString, jobname::AbstractString, datadir::AbstractString=joinpath(workdir, "data"), SimID::Int=0, MPIrank::Int=0, resuming::Bool=true)
    jobname_str = (SimID==0 ? string(jobname) : join([jobname, "$(SimID)"], "_"))
    maindir = joinpath(datadir, jobname_str)

    if resuming # default setting
        resuming = isdir(maindir)
    else # if you insist on replacing the old data
        if isdir(maindir)
            @warn "The directory $maindir already exists and it will be replaced."
        end
    end

    return SimulationInfo(string(workdir), string(datadir), jobname_str, maindir, MPIrank, SimID, resuming)
end

"""
    SimulationInfo(prefix::AbstractString, paras_info::Vector{<:Tuple{String,Any,Int}}; workdir::AbstractString, kwargs...)

Create a SimulationInfo object by generating jobname from prefix and parameters.

# Arguments
- `prefix::AbstractString`: The prefix name for the job (e.g., "honeycomb")
- `paras_info::Vector{<:Tuple{String,Any,Int}}`: Vector of parameter tuples (key, value, digits)
- `workdir::AbstractString`: The path to the project directory
- `kwargs...`: Additional keyword arguments passed to the main constructor

# Returns
- `SimulationInfo`: The constructed simulation info object

# Examples
```julia
# Mixed parameter types
params = [("L", 8, 0), ("U", 4.0, 1), ("enable_afm", true, 0)]
sim_info = SimulationInfo("honeycomb", params, workdir="/path/to/project")
```
"""
function SimulationInfo(prefix::AbstractString, paras_info::Vector{<:Tuple{String,Any,Int}}; workdir::AbstractString, kwargs...)
    jobname = generate_jobname(prefix, paras_info)
    return SimulationInfo(workdir=workdir, jobname=jobname; kwargs...)
end

"""
    Base.show(io::IO, ::MIME"text/plain", sim_info::SimulationInfo)

Display SimulationInfo object in TOML format.

# Arguments
- `io::IO`: Output stream
- `::MIME"text/plain"`: MIME type for plain text display
- `sim_info::SimulationInfo`: The SimulationInfo object to display
"""
function Base.show(io::IO, ::MIME"text/plain", sim_info::SimulationInfo)

    @printf io "[simulation_info]\n\n"
    @printf io "jobname  = \"%s\"\n" sim_info.jobname
    @printf io "SimID    = %d\n" sim_info.SimID
    @printf io "MPIrank  = %d\n" sim_info.MPIrank
    @printf io "julia_version     = \"%s\"" VERSION

    return nothing
end

"""
    save_simulation_info(sim_info::SimulationInfo, additional_info=nothing)

Save simulation information to a TOML file in the main directory.

# Arguments
- `sim_info::SimulationInfo`: The simulation info object to save
- `additional_info`: Optional additional information to include in the TOML file (default: nothing)

# Examples
```julia
save_simulation_info(sim_info)
save_simulation_info(sim_info, Dict("custom_param" => "value"))
```
"""
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