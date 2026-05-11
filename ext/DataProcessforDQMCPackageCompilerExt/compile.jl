using DataProcessforDQMC: DataProcessforDQMC, Algorithm
using PackageCompiler: PackageCompiler

"""
    compile(::Algorithm{:PackageCompiler}; dir, filename, additional_packages=Symbol[])

构建包含 `DataProcessforDQMC` 的预编译系统镜像。

# Keyword Arguments
- `dir`, `filename`: sysimage 的输出目录与文件名。
- `additional_packages`: 额外打包进 sysimage 的包名列表（例如 `[:OhMyREPL, :Revise]`）。
  这些包须已经安装在当前 active project 中。把启动时 `using` 的常用包打进 sysimage
  可以避免 `--sysimage` 启动后再触发 precompile —— 个别包（如 OhMyREPL v0.5.32）
  在带 sysimage 的 precompile 子进程里会因 Pkg `REPLExt` extension 不激活而失败，
  baked-in 是最干净的规避方式。

# Example
```julia
DataProcessforDQMC.compile(additional_packages=[:OhMyREPL])
```
"""
function DataProcessforDQMC.compile(
  ::Algorithm{:PackageCompiler};
  dir::AbstractString=DataProcessforDQMC.default_compile_dir(),
  filename::AbstractString=DataProcessforDQMC.default_compile_filename(),
  additional_packages=Symbol[],
)
  if !isdir(dir)
    println("""The directory "$dir" doesn't exist yet, creating it now.""")
    println()
    mkdir(dir)
  end
  path = joinpath(dir, filename)
  println(
    """Creating the system image "$path" containing the compiled version of DataProcessforDQMC. This may take a few minutes.""",
  )
  packages = isempty(additional_packages) ?
             :DataProcessforDQMC :
             Symbol[:DataProcessforDQMC; collect(Symbol, additional_packages)]
  PackageCompiler.create_sysimage(
    packages;
    sysimage_path=path,
    precompile_execution_file=joinpath(@__DIR__, "precompile_dataprocessfordqmc.jl"),
  )
  println(DataProcessforDQMC.compile_note(; dir, filename))
  return path
end