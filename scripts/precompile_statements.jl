# scripts/precompile_statements.jl
# 这个文件模仿 run_scnet.jl 的核心调用，以触发关键函数的预编译。

using DataProcessforDQMC.JobManage

# 将预编译调用包裹在 try-catch 块中是很好的实践。
# 这样即使因为文件不存在等原因而出错，编译本身也已经触发了。
try
    # --- 模仿 run_scnet.jl 和 cal_para.jl 中的设置 ---

    # 1. 定义一些典型的参数
    code = "proj_fft_honeycomb_exact"
    beta = 4.0
    u = -3.50
    L = 3
    dtau = 0.05
    workdir = "."
    datadir = "."

    # 2. 构建与脚本中类型一致的参数数组
    job_para_info = [("b", beta, 3), ("U", u, 2), ("L", L, 0), ("dtau", dtau, 2)]

    # 3. 【核心】调用关键函数来触发其编译
    # 这个调用会编译 SimulationInfo 构造函数，使其在未来能够被快速调用。
    SimInfo = SimulationInfo(code, job_para_info; workdir=workdir, datadir=datadir)

catch e
    # 这里的失败通常是可以接受的，我们的目标是触发编译，而不是成功运行。
    @warn "Precompilation statements failed (this is often okay): $e"
end 