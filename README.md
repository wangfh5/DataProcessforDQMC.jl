# DataProcessforDQMC

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wangfh5.github.io/DataProcessforDQMC.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wangfh5.github.io/DataProcessforDQMC.jl/dev/)
[![Build Status](https://github.com/wangfh5/DataProcessforDQMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wangfh5/DataProcessforDQMC.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/wangfh5/DataProcessforDQMC.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/wangfh5/DataProcessforDQMC.jl)

## Intention of the package

Provide tools for processing the data (various `.bin` files) output by DQMC codes. 

- `DataProcessforDQMC.jl` does not care about the model parameters at most time. Instead, it accepts the location directory(s) of the data you would like to format or read. 