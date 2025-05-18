#= 
参数分析脚本 (Parameter Analysis)

此脚本用于在包含多个参数文件夹的目录中运行分析，可以:
1. 对单个参数进行分析，返回数组
2. 对多个参数进行分析，返回DataFrame

使用示例:
1. 单参数分析: 
   julia -e 'using DataProcessforDQMC; analyze_parameter(".", "U", 4.0)'

2. 多参数分析:
   julia -e 'using DataProcessforDQMC; analyze_multiple_parameters(".", "U", [2.0, 3.0, 4.0])'
=#
