# [预编译系统镜像指南](@id precompilation_guide)

为了获得最佳的启动速度和性能，DataProcessforDQMC支持预编译系统镜像。预编译后的Julia启动时间从数秒缩短到毫秒级。

## 🚀 前置准备（首次使用）

如果您还没有安装Julia，请按照以下步骤操作：

### 1. 安装Julia

使用juliaup安装器（官方推荐）：

```bash
curl -fsSL https://install.julialang.org | sh
```

安装完成后，重新加载shell配置：
```bash
source ~/.bashrc  # 或 source ~/.zshrc
```

验证安装：
```bash
julia --version
```

### 2. 安装PackageCompiler

启动Julia REPL并安装PackageCompiler包：

```bash
julia
```

在Julia REPL中执行：
```julia
using Pkg
Pkg.add("PackageCompiler")
```

### 3. 安装DataProcessforDQMC包到全局环境

在Julia REPL中安装本地包到全局环境：

```bash
julia
```

在Julia REPL中执行：

```julia
using Pkg

# 方式1：仅使用（推荐）
Pkg.add(path="/path/to/DataProcessforDQMC")

# 方式2：开发模式（如需修改源代码）
Pkg.develop(path="/path/to/DataProcessforDQMC")
```

!!! tip "安装方式选择"
    - 使用 `Pkg.add()` 如果您只需要使用这个包
    - 使用 `Pkg.develop()` 如果您需要修改包的源代码

## 📦 构建系统镜像

启动 Julia 并执行：

```julia
using DataProcessforDQMC, PackageCompiler
DataProcessforDQMC.compile()
```

这将在 `~/.julia/sysimages/` 目录下生成 `sys_dataprocessfordqmc.so` 系统镜像文件（约380MB），并自动显示使用说明。

!!! note "编译时间"
    首次编译可能需要1-2分钟，请耐心等待。

### 把额外的常用包打进 sysimage

`compile()` 接受 `additional_packages` 关键字参数，把启动时常用的其它包（如 `OhMyREPL`、`Revise`、`Plots` 等）一并 baked 进 sysimage，避免 `jd` 启动时再触发 precompile：

```julia
DataProcessforDQMC.compile(additional_packages=[:OhMyREPL])
```

列表中的包必须已经安装在当前 active project 中（通常是 `@v1.x` 默认环境）。

!!! tip "OhMyREPL 强烈建议加入"
    Julia 1.11 + OhMyREPL v0.5.32 组合下，若用户 `startup.jl` 在 `atreplinit` 里 `using OhMyREPL`，启动 `jd` 时会触发一次 precompile 子进程；该子进程继承 sysimage 后，Pkg 的 `REPLExt` extension 不会自动激活，`OhMyREPL/src/BracketInserter.jl` 会以 `type Nothing has no field promptf` 报错。把 `:OhMyREPL` 加入 `additional_packages` 让它在 sysimage 中直接 baked 加载，可彻底规避该问题。

## ⚙️ 配置Julia启动器

将以下配置添加到您的 `~/.bashrc` 或 `~/.zshrc` 文件中：

```bash
# DataProcessforDQMC 预编译镜像启动器
function jd() {
    local sysimage="$HOME/.julia/sysimages/sys_dataprocessfordqmc.so"

    if [ ! -f "$sysimage" ]; then
        echo "Error: system image file not found: $sysimage"
        return 1
    fi

    if [ $# -eq 0 ]; then
        # 无参数：启动交互式会话并自动加载包
        echo "Starting DataProcessforDQMC interactive session..."
        julia --sysimage "$sysimage" -e 'using DataProcessforDQMC' -i
    else
        # 有参数：运行脚本或传递其他参数
        julia --sysimage "$sysimage" "$@"
    fi
}

# 开发版本（不使用系统镜像，实时加载代码修改）
alias jdd="julia -e 'using DataProcessforDQMC' -i"
```

然后重新加载配置：
```bash
source ~/.bashrc  # 或 source ~/.zshrc
```

## 🎯 使用方式

### 交互式会话
```bash
jd  # 启动预编译的DataProcessforDQMC交互会话
```

### 运行脚本
```bash
jd your_analysis_script.jl  # 使用预编译镜像运行脚本
```

### 开发调试
```bash
jdd  # 使用开发版本（实时加载代码修改）
```

## 🔄 何时需要重新编译

以下情况需要重新编译系统镜像：

1. **代码修改**：修改了DataProcessforDQMC包的源代码
2. **依赖更新**：更新了包的依赖项
3. **Julia版本升级**：升级了Julia版本
4. **函数签名变更**：修改了公共API函数的参数

重新编译命令：
```julia
using DataProcessforDQMC, PackageCompiler
DataProcessforDQMC.compile()
```

!!! warning "代码修改后"
    每次修改源代码后都需要重新编译系统镜像，否则修改不会生效。开发调试时建议使用 `jdd` 命令，它会实时加载代码修改。

## 💡 性能对比

| 启动方式 | 首次加载时间 | 后续启动时间 | 适用场景 |
|---------|-------------|-------------|----------|
| `jd` (预编译) | ~200ms | ~100ms | 生产分析、批量处理 |
| `jdd` (开发版) | ~3-5s | ~2-3s | 代码开发、调试 |

## 📋 基本使用示例

### 多参数分析
```julia
# 启动预编译版本
jd

# 分析反铁磁相关比率
base_dir = "/path/to/simulation/data"
results = analyze_AFM_correlation_ratio_multi_parameter(
    base_dir,
    shift_point=(0.25, 0.0),
    filter_options=Dict(
        :prefix => "proj_bt_honeycomb_exact",
        :L => 6,
        :dtau => 0.1
    )
)
```

### 目录扫描
```julia
# 扫描符合条件的参数目录
dirs = scan_parameter_directories(
    base_dir,
    filter_options=Dict(
        :prefix => "honeycomb_model",
        :U => [3.0, 4.0, 5.0],
        :L => 12
    )
)
```

## 🛠️ 故障排除

### 系统镜像损坏
如果遇到奇怪的错误，尝试重新编译：
```bash
rm ~/.julia/sysimages/sys_dataprocessfordqmc.so
```
然后重新运行：
```julia
using DataProcessforDQMC, PackageCompiler
DataProcessforDQMC.compile()
```

### 路径问题
确认系统镜像文件已正确生成：
```bash
ls -la ~/.julia/sysimages/sys_dataprocessfordqmc.so
```

### PackageCompiler未安装
如果编译失败，可能需要安装PackageCompiler：
```julia
using Pkg
Pkg.add("PackageCompiler")
```

### `jd` 启动时反复触发 OhMyREPL precompile 并失败
现象：每次 `jd` 启动都出现 `Precompiling OhMyREPL...` → `type Nothing has no field promptf` 错误。

原因：Julia 1.11 在带 sysimage 的 precompile 子进程里不会激活 Pkg 的 `REPLExt` extension，OhMyREPL 顶层依赖该 extension 来取 `promptf`。

解决：把 OhMyREPL baked 进 sysimage：
```julia
using DataProcessforDQMC, PackageCompiler
DataProcessforDQMC.compile(additional_packages=[:OhMyREPL])
```

!!! tip "建议"
    建议将预编译版本 `jd` 用于日常分析工作，将开发版本 `jdd` 用于代码开发和调试。
