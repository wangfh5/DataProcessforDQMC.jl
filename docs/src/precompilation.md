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

在项目根目录执行：

```bash
cd /path/to/DataProcessforDQMC
julia --project=. scripts/build_sysimage.jl
```

这将生成 `DataProcessforDQMC_sys.so` 系统镜像文件（约380MB）。

!!! note "编译时间"
    首次编译可能需要1-2分钟，请耐心等待。

## ⚙️ 配置Julia启动器

将以下配置添加到您的 `~/.bashrc` 文件中：

```bash
# DataProcessforDQMC 预编译镜像启动器
function julia_dqmc() {
    local sysimage_path="/path/to/DataProcessforDQMC/DataProcessforDQMC_sys.so"
    
    if [ ! -f "$sysimage_path" ]; then
        echo "Error: system image file not found: $sysimage_path"
        return 1
    fi
    
    export JULIA_DQMC_SESSION=true
    
    if [ $# -eq 0 ]; then
        echo "Starting DataProcessforDQMC interactive session..."
        julia -J"$sysimage_path" -e "using DataProcessforDQMC" -i
    else
        julia -J"$sysimage_path" "$@"
    fi
    
    unset JULIA_DQMC_SESSION
}

# 快捷命令
alias jd='julia_dqmc'        # 预编译版本
alias jdd="julia -e \"using DataProcessforDQMC\" -i "  # 开发版本
```

然后重新加载配置：
```bash
source ~/.bashrc
```

!!! tip "路径配置"
    请将 `/path/to/DataProcessforDQMC` 替换为实际的项目路径。

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
```bash
cd /path/to/DataProcessforDQMC
julia --project=. scripts/build_sysimage.jl
```

!!! warning "代码修改后"
    每次修改源代码后都需要重新编译系统镜像，否则修改不会生效。

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
rm DataProcessforDQMC_sys.so
julia --project=. scripts/build_sysimage.jl
```

### 路径问题
确保 `.bashrc` 中的路径指向正确的系统镜像文件：
```bash
ls -la /path/to/DataProcessforDQMC/DataProcessforDQMC_sys.so
```

### 权限问题
确保系统镜像文件有执行权限：
```bash
chmod +x DataProcessforDQMC_sys.so
```

### PackageCompiler问题
如果编译失败，可能需要安装或更新PackageCompiler：
```julia
using Pkg
Pkg.add("PackageCompiler")
```

!!! tip "建议"
    建议将预编译版本 `jd` 用于日常分析工作，将开发版本 `jdd` 用于代码开发和调试。
