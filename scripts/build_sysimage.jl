using PackageCompiler

# 将你希望包含在系统镜像中的包列在这里
# 核心包是 DataProcessforDQMC
const packages_to_compile = [:DataProcessforDQMC]

# 预编译语句文件的路径
const precompile_statements_file = joinpath(@__DIR__, "precompile_statements.jl")

# 输出的系统镜像文件名
# 注意：在Linux/macOS上后缀是 .so, 在Windows上是 .dll
const sysimage_path = joinpath(@__DIR__, "..", "DataProcessforDQMC_sys.so")

println("Starting system image build...")

create_sysimage(
    packages_to_compile,
    sysimage_path=sysimage_path,
    precompile_execution_file=precompile_statements_file
)

println("System image build complete!")
println("Image saved to: ", sysimage_path) 