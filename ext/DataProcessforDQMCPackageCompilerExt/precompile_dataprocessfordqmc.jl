# Precompile statements for DataProcessforDQMC
# This file executes typical usage patterns to trigger compilation

using DataProcessforDQMC.JobManage

# Wrap precompile calls in try-catch block
# Even if execution fails (e.g., due to missing files), compilation is triggered
try
    # --- Simulate typical usage from run_scnet.jl and cal_para.jl ---

    # 1. Define typical parameters
    code = "proj_fft_honeycomb_exact"
    beta = 4.0
    u = -3.50
    L = 3
    dtau = 0.05
    workdir = "."
    datadir = "."

    # 2. Build parameter array with consistent types
    job_para_info = [("b", beta, 3), ("U", u, 2), ("L", L, 0), ("dtau", dtau, 2)]

    # 3. [KEY] Call critical functions to trigger their compilation
    # This compiles the SimulationInfo constructor for faster future calls
    SimInfo = SimulationInfo(code, job_para_info; workdir=workdir, datadir=datadir)

catch e
    # Failure here is often acceptable - our goal is to trigger compilation, not successful execution
    @warn "Precompilation statements failed (this is often okay): $e"
end