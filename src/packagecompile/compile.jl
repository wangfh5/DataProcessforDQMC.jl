# Algorithm dispatch type for compile backends
struct Algorithm{Name} end
Algorithm(s::String) = Algorithm{Symbol(s)}()
macro Algorithm_str(s)
  return :(Algorithm{$(Expr(:quote, Symbol(s)))}())
end

default_compile_dir() = joinpath(homedir(), ".julia", "sysimages")

default_compile_filename() = "sys_dataprocessfordqmc.so"

default_compile_path() = joinpath(default_compile_dir(), default_compile_filename())

function compile_note(; dir=default_compile_dir(), filename=default_compile_filename())
  path = joinpath(dir, filename)
  return """
  ✓ 系统镜像编译成功！

  将以下配置添加到您的 ~/.bashrc 或 ~/.zshrc 文件中：

    # DataProcessforDQMC 预编译镜像启动器
    function jd() {
        local sysimage="\$HOME/.julia/sysimages/sys_dataprocessfordqmc.so"

        if [ ! -f "\$sysimage" ]; then
            echo "Error: system image file not found: \$sysimage"
            return 1
        fi

        if [ \$# -eq 0 ]; then
            # 无参数：启动交互式会话并自动加载包
            echo "Starting DataProcessforDQMC interactive session..."
            julia --sysimage "\$sysimage" -e 'using DataProcessforDQMC' -i
        else
            # 有参数：运行脚本或传递其他参数
            julia --sysimage "\$sysimage" "\$@"
        fi
    }

    # 开发版本（不使用系统镜像，实时加载代码修改）
    alias jdd="julia -e 'using DataProcessforDQMC' -i"

  然后运行 source ~/.bashrc 即可使用：
    jd              # 交互式会话
    jd script.jl    # 运行脚本
  """
end

function compile(; backend=Algorithm"PackageCompiler", kwargs...)
  return compile(backend; kwargs...)
end