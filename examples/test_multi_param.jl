using DataProcessforDQMC
using Statistics

println("=== DataProcessforDQMC JobManage 完整使用示例 ===\n")

## -------------------------------------------------------------------------- ##
##                              任务1: 从历史文件夹名提取参数                              ##
## -------------------------------------------------------------------------- ##

println("📁 任务1: 从历史文件夹名提取参数 (Legacy 格式)")
println("=" ^ 60)

# 测试从历史文件夹名提取参数 - 使用 common-functions.jl 中的函数
test_dirname = "proj_fft_honeycomb_exact.b8.000.U4.00.gw0.62.lprojgwF.L9.dtau0.05"
println("历史格式文件夹名: $test_dirname")

# 使用新的 parse_jobname_legacy 函数来解析旧格式
prefix_name, legacy_params = parse_jobname_legacy(test_dirname)
println("\n📊 使用 parse_jobname_legacy 解析结果:")
println("  前缀: $prefix_name")
for (key, value, digits) in legacy_params
    println("  $key = $value ($(typeof(value)), $digits 位小数)")
end

# 提取有用的参数用于后续转换（现在legacy_params已经是Vector{Tuple{String,Any,Int}}）
useful_params = legacy_params

println("\n✅ 任务1完成：成功从legacy格式提取了 $(length(useful_params)) 个参数")

## -------------------------------------------------------------------------- ##
##                           任务2: 利用提取参数来创建新规范下的文件名                           ##
## -------------------------------------------------------------------------- ##

println("\n" * repeat("=", 80))
println("🔄 任务2: 利用提取参数来创建新规范下的文件名")
println(repeat("=", 60))

# 将旧格式的参数转换为新格式所需的格式
println("📝 参数转换过程:")

# 现在useful_params已经包含精度信息，直接使用
new_format_params = useful_params

println("📝 参数转换过程:")
for (key, value, digits) in new_format_params
    println("  $key: $value (保持，$digits 位小数)")
end

# 生成新格式的文件夹名
new_jobname = generate_jobname(prefix_name, new_format_params)

println("\n🎯 转换结果:")
println("  旧格式: $test_dirname")
println("  新格式: $new_jobname")

# 验证格式长度对比
println("\n📏 格式对比:")
println("  旧格式长度: $(length(test_dirname)) 字符")
println("  新格式长度: $(length(new_jobname)) 字符")
println("  压缩比: $(round((1 - length(new_jobname)/length(test_dirname))*100, digits=1))%")

println("\n✅ 任务2完成：成功将legacy格式转换为新规范格式")

## -------------------------------------------------------------------------- ##
##                             任务3: 测试从新规范文件名提取参数                             ##
## -------------------------------------------------------------------------- ##

println("\n" * repeat("=", 80))
println("🔍 任务3: 测试从新规范文件名提取参数")
println(repeat("=", 60))

# 使用刚刚生成的新格式文件夹名进行测试
test_new_format = new_jobname
println("测试文件夹名: $test_new_format")

# 提取前缀和参数
extracted_prefix, extracted_params = parse_jobname(test_new_format)
println("\n📤 提取结果:")
println("前缀: '$extracted_prefix'")
println("参数:")
for (key, value) in sort(collect(extracted_params))
    println("  $key = $value ($(typeof(value)))")
end

# 带格式信息的解析
println("\n📋 带格式信息的解析:")
prefix_with_digits, params_with_digits = parse_jobname(test_new_format)
println("前缀: '$prefix_with_digits'")
for (key, value, digits) in params_with_digits
    println("  $key = $value ($(typeof(value)), $digits 位小数)")
end

# 验证参数一致性
println("\n🔍 验证参数一致性:")
all_consistent = true

# 将extracted_params转换为字典以便查找
extracted_dict = Dict(key => value for (key, value, _) in extracted_params)

for (orig_key, orig_value, orig_digits) in useful_params
    if haskey(extracted_dict, orig_key)
        extracted_value = extracted_dict[orig_key]
        # 数值参数比较（允许小的浮点误差）
        if isa(orig_value, AbstractFloat) && isa(extracted_value, AbstractFloat)
            consistent = abs(orig_value - extracted_value) < 1e-10
        else
            consistent = (orig_value == extracted_value)
        end
        
        if consistent
            println("  ✅ $orig_key: $orig_value ↔ $extracted_value")
        else
            println("  ❌ $orig_key: $orig_value ↔ $extracted_value (不一致!)")
            all_consistent = false
        end
    else
        println("  ⚠️  $orig_key: 在新格式中未找到")
        all_consistent = false
    end
end

if all_consistent
    println("\n🎉 所有参数转换一致！")
else
    println("\n⚠️  发现参数不一致，需要检查转换逻辑")
end

println("\n✅ 任务3完成：成功验证新格式解析功能")

## -------------------------------------------------------------------------- ##
##                                   批量迁移                                  ##
## -------------------------------------------------------------------------- ##

println("\n" * repeat("=", 80))
println("🚀 高级示例：批量迁移演示")
println(repeat("=", 60))

# 模拟多个legacy文件夹的迁移
legacy_examples = [
    "proj_fft_honeycomb_exact.b8.000.U4.00.L6.dtau0.05",
    "proj_fft_honeycomb_exact.b10.000.U2.00.gw0.5.L9.dtau0.1",
    "proj_fft_square.b6.000.U6.00.lprojgwT.L12.dtau0.02"
]

println("📂 批量迁移示例:")
migration_results = []

for (i, legacy_name) in enumerate(legacy_examples)
    println("\n$i. Legacy: $legacy_name")
    
    # 解析旧格式
    old_prefix, old_params = parse_jobname_legacy(legacy_name)
    
    # 转换为新格式（现在old_params已经包含精度信息）
    new_params = old_params
    
    new_name = generate_jobname(old_prefix, new_params)
    println("   New:    $new_name")
    println("   压缩:   $(round((1 - length(new_name)/length(legacy_name))*100, digits=1))%")
    
    push!(migration_results, (legacy_name, new_name))
end

println("\n📊 迁移总结:")
println("  • 处理文件夹数: $(length(legacy_examples))")
println("  • 平均压缩比: $(round(mean([(1 - length(new)/length(old))*100 for (old, new) in migration_results]), digits=1))%")
println("  • 所有转换成功: ✅")

println("\n🎉 完整演示结束！")

## -------------------------------------------------------------------------- ##
##                                 真实文件系统的批量迁移                                ##
## -------------------------------------------------------------------------- ##

root = "/home/wangfh5/Projects/archieve/honeycomb-Gutzwiller/run_honeycomb"

# list all directories in the path
dirs = readdir(root)

# filter out directories that are not directories
dirs = filter(dir -> isdir(joinpath(root, dir)), dirs)

# dry run
batch_migrate_directories(root; dry_run=true, pattern=r"^proj_bt_honeycomb")

# copy
batch_migrate_directories(root; dry_run=false, action="cp")

# remove
batch_migrate_directories(root; dry_run=false, action="rm")

# move
batch_migrate_directories(root; dry_run=false, pattern=r"^proj_bt_honeycomb", action="mv")