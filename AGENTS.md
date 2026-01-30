---
description: 
alwaysApply: true
---

# Repository Guidelines

## Project Structure & Module Organization
- `src/` core package code. Submodules: `JobManage/`, `data-processing/`, `single-parameter-analysis/`, `multiple-parameter-analysis/`, `packagecompile/`.
- `test/` unit tests (`runtests.jl`). Add feature tests as separate files and include them from `runtests.jl`.
- `docs/` Documenter.jl site (`docs/make.jl`).
- `examples/` runnable scripts (e.g., `examples/test_multi_param.jl`).
- `ext/` Package extension for optional `PackageCompiler` integration.

## Build, Test, and Development Commands
- Setup env: `julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'`
- Run tests (with coverage): `julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'`
- Run example: `julia --project=. examples/test_multi_param.jl`
- Build docs locally: `julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'`
- Optional sysimage (faster startup): See detailed instructions in `docs/src/precompilation.md` (requires global environment setup with PackageCompiler).

## Coding Style & Naming Conventions
- Julia style, 4‑space indentation; keep lines readable (<100 chars); no trailing whitespace.
- Names: Modules/Types `TitleCase`, functions `snake_case`, constants `UPPER_CASE`.
- Follow existing file naming (hyphenated, e.g., `structure-factor.jl`).
- Add docstrings for exported APIs; keep exports centralized in `src/DataProcessforDQMC.jl`.

## Testing Guidelines
- Use `Test` stdlib. Organize tests by feature with `@testset` blocks.
- Create files like `test_structure_factor.jl` and include from `test/runtests.jl`.
- Prefer deterministic tests; seed randomness when unavoidable.
- Coverage is reported in CI (Codecov). Don’t decrease coverage; add tests for new or changed public APIs.

## Data & Repository Hygiene
- Do not commit large DQMC outputs or generated binaries; keep datasets local or in releases.
- Place small, illustrative inputs under `examples/` if needed; prefer reproducible scripts over raw data dumps.
