# 文档部署指南

本文档说明如何将 DataProcessforDQMC.jl 的文档部署到 GitHub Pages。

## 自动部署流程概览

```
代码推送到 main 分支
    ↓
GitHub Actions 自动触发
    ↓
CI 运行 docs/make.jl
    ↓
Documenter.jl 构建 HTML
    ↓
自动推送到 gh-pages 分支
    ↓
GitHub Pages 自动发布
    ↓
文档可在 https://wangfh5.github.io/DataProcessforDQMC.jl 访问
```

## 一次性设置步骤

### 1. 启用 GitHub Pages（仓库设置）

1. 访问 https://github.com/wangfh5/DataProcessforDQMC.jl/settings/pages
2. 在 **Source** 部分：
   - Branch: 选择 `gh-pages`
   - Folder: 选择 `/ (root)`
3. 点击 **Save**
4. 等待几分钟，页面会显示：
   ```
   ✅ Your site is live at https://wangfh5.github.io/DataProcessforDQMC.jl/
   ```

### 2. 生成部署密钥

在本地运行（**只需运行一次**）：

```bash
cd /Users/ssqc/Projects/DataProcessforDQMC.jl

# 安装 DocumenterTools（如果还没有）
julia -e 'using Pkg; Pkg.add("DocumenterTools")'

# 生成密钥
julia docs/generate_deploy_key.jl
```

**重要提示：**
- 这个脚本会生成一对 SSH 密钥
- **公钥**会自动添加到 GitHub 仓库的 Deploy Keys
- **私钥**需要你手动复制到 GitHub Secrets（见下一步）

### 3. 添加私钥到 GitHub Secrets

1. 复制上一步输出的私钥（完整内容，包括 `-----BEGIN OPENSSH PRIVATE KEY-----` 等）
2. 访问 https://github.com/wangfh5/DataProcessforDQMC.jl/settings/secrets/actions
3. 点击 **New repository secret**
4. 填写：
   - Name: `DOCUMENTER_KEY`
   - Value: [粘贴私钥]
5. 点击 **Add secret**

### 4. 验证 CI 配置

确认文件 `.github/workflows/CI.yml` 中包含 `docs` job（已配置）。

## 使用方法

### 自动部署（推荐）

每次推送到 `main` 分支时，GitHub Actions 会自动：
1. 构建文档
2. 部署到 `gh-pages` 分支
3. 更新 https://wangfh5.github.io/DataProcessforDQMC.jl/

**查看构建状态：**
- 访问 https://github.com/wangfh5/DataProcessforDQMC.jl/actions
- 查看 "CI" workflow 中的 "Documentation" job 运行状态

### 本地预览（开发时）

**方法 1：普通构建**
```bash
julia docs/make.jl
open docs/build/index.html  # 用 file:// 协议打开（导航链接可能有问题）
```

**方法 2：带预览服务器（推荐）**
```bash
DOCS_PREVIEW=true julia docs/make.jl
# 自动启动 http://localhost:8000
# 浏览器访问 http://localhost:8000
# 按 Ctrl+C 停止
```

**方法 3：手动启动服务器**
```bash
# 先构建
julia docs/make.jl

# 再启动服务器
cd docs/build
python3 -m http.server 8000

# 浏览器访问 http://localhost:8000
```

## 故障排查

### 问题 1: CI 构建失败

**查看日志：**
1. 访问 https://github.com/wangfh5/DataProcessforDQMC.jl/actions
2. 点击失败的 workflow run
3. 查看详细错误信息

**常见原因：**
- `DOCUMENTER_KEY` 未设置或设置错误
- 文档源文件中有语法错误
- 交叉引用 `@ref` 失效

### 问题 2: 部署成功但页面显示 404

**检查 GitHub Pages 设置：**
1. 访问 https://github.com/wangfh5/DataProcessforDQMC.jl/settings/pages
2. 确认 Source 设置为 `gh-pages` 分支
3. 等待 2-5 分钟让 GitHub 重新部署

### 问题 3: 本地预览导航链接不工作

**原因：** 使用 `file://` 协议打开时，浏览器不会自动加载目录下的 `index.html`

**解决：** 使用 HTTP 服务器预览（见上文"本地预览"）

## 文档结构

```
docs/
├── make.jl                    # 构建脚本
├── Project.toml               # docs 环境依赖
├── generate_deploy_key.jl     # 密钥生成脚本（一次性使用）
├── DEPLOYMENT.md              # 本文档
└── src/
    ├── index.md               # 首页
    ├── precompilation.md      # 预编译指南
    ├── multi_parameter_analysis.md  # 多参数分析
    ├── bin_analysis.md        # Bin 分析
    └── api.md                 # API 参考

.github/workflows/
└── CI.yml                     # CI 配置（包含 docs job）
```

## 相关链接

- **文档站点**: https://wangfh5.github.io/DataProcessforDQMC.jl/
- **GitHub Actions**: https://github.com/wangfh5/DataProcessforDQMC.jl/actions
- **Documenter.jl 文档**: https://documenter.juliadocs.org/stable/

