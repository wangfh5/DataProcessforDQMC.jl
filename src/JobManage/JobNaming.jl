"""
    generate_jobname(prefix::AbstractString, paras_info::Vector{Tuple{String,Any,Int}})

The definition of a jobname: 
1. start with the prefix name `prefix`, e.g., "free_fermion";
2. followed by all the essential parameters, e.g., "L=8_t=1.0_mu=0.500".
Each element in `paras_info` is a tuple, which contains the `(key, value, digits)` of a parameter. 

The design of the jobname:
- `=` is chosen as the delimiter between the key and value of each parameter (learned from ALF).
- `_` is chosen as the delimiter between parameter strings (e.g., `L=8` and `t=1.0` are separated by `_`), 
since `-` and `.` could present in the parameter values.

Some rules of the `key` of parameter and the prefix name `prefix`:
1. Do NOT contain the delimiters `_` and `=` in the `keys`. 
2. Do NOT contain the delimeter `=` in the `prefix`.
3. We recommend that no special unicode characters are used, e.g., `λ` should be replaced by `lambda`. 
Since it is hard to type the unicode characters in the terminal.

Supports Bool parameters: Bool true becomes "T", Bool false becomes "F".
For Real values, digits specifies decimal places.

# Examples
```julia
# Mixed parameter types
params = [("L", 8, 0), ("U", 4.0, 1), ("enable_afm", true, 0)]
jobname = generate_jobname("honeycomb", params)
# Result: "honeycomb_L=8_U=4.0_enable_afm=T"
```
"""
function generate_jobname(prefix::AbstractString, paras_info::Vector{<:Tuple{String,Any,Int}})
    jobname = prefix
    for para_info ∈ paras_info
        key, value, digits = para_info
        
        if value isa Bool
            # Convert Bool to T/F following the design philosophy
            value_str = value ? "T" : "F"
        elseif value isa Real
            # Format Real numbers with specified digits
            fmt = Printf.Format("%.$(digits)f")
            value_str = Printf.format(fmt, value)
        else
            # For other types, convert to string directly
            value_str = string(value)
        end
        
        parastr = join([key, value_str], "=")
        jobname = join([jobname, parastr], "_")
    end
    return jobname
end

"""
    parse_jobname(jobname::AbstractString)

Parse a jobname and extract both prefix and parameters.
This is the recommended function for parsing jobnames.

# Arguments
- `jobname::AbstractString`: The jobname to parse

# Returns
- `(prefix::String, params::Vector{Tuple{String,Any,Int}})`: prefix and parameters with digits info

# Examples
```julia
prefix, params = parse_jobname("honeycomb_L=8_U=4.0_enable_afm=T")
# prefix = "honeycomb"
# params = [("L", 8, 0), ("U", 4.0, 1), ("enable_afm", true, 0)]
```
"""
function parse_jobname(jobname::AbstractString)
    # 1. 找到第一个等号的位置
    first_eq_pos = findfirst('=', jobname)
    
    if first_eq_pos === nothing
        # 没有等号，整个字符串都是前缀
        return jobname, Tuple{String,Any,Int}[]
    end
    
    # 2. 提取第一个等号之前的所有内容
    raw_prefix = jobname[1:first_eq_pos-1]
    
    # 3. 找到最后一个下划线，去掉该下划线后的参数名
    last_underscore_pos = findlast('_', raw_prefix)
    
    if last_underscore_pos === nothing
        # 没有下划线，这种情况很少见，可能是直接以参数开头
        prefix = ""
    else
        prefix = jobname[1:last_underscore_pos-1]
    end
    
    # 4. 用正则表达式匹配所有参数
    jobname_without_prefix = jobname[length(prefix)+1:end]
    all_matches = collect(eachmatch(r"_([a-zA-Z_]+)=([TF]|[\d\.\-]+)", jobname_without_prefix))
    
    # 如果没有匹配到参数
    if isempty(all_matches)
        return prefix, Tuple{String,Any,Int}[]
    end
    
    # 5. 解析参数
    paras_info = Vector{Tuple{String,Any,Int}}()
    for match_obj in all_matches
        key, value_str = match_obj.captures
        value, digits = _parse_parameter_value(value_str)
        push!(paras_info, (key, value, digits))
    end
    return prefix, paras_info
end

"""
    parse_jobname_legacy(jobname::AbstractString)

从目录名中提取参数值，包括前缀类型、温度参数b、相互作用U、晶格尺寸L和时间步长dtau。
支持legacy格式的参数解析，如 "proj_fft_honeycomb_exact.b8.000.U4.00.L9.dtau0.05"

# 参数
- `jobname::AbstractString`: 目录名，例如 "proj_fft_honeycomb_exact.b8.000.U4.00.L9.dtau0.05"

# 返回值
- `(prefix::String, params::Vector{Tuple{String,Any,Int}})`: 前缀和参数列表（保持原始顺序）

# 示例
```julia
prefix, params = parse_jobname_legacy("proj_fft_honeycomb_exact.b8.000.U4.00.L9.dtau0.05")
# prefix = "proj_fft_honeycomb_exact"
# params = [("b", 8.0, 3), ("U", 4.0, 2), ("L", 9, 0), ("dtau", 0.05, 2)]
```
"""
function parse_jobname_legacy(jobname::AbstractString)
    # 使用正则表达式提取前缀，支持更多格式
    prefix_patterns = [
        r"^(proj_fft_honeycomb(?:_[a-zA-Z]+)?)",  # honeycomb格式
        r"^(proj_fft_square(?:_[a-zA-Z]+)?)",     # square格式
        r"^(proj_fft_chain(?:_[a-zA-Z]+)?)",      # chain格式
        r"^(proj_fft_cubic(?:_[a-zA-Z]+)?)",      # cubic格式
        r"^([a-zA-Z_]+)(?=\.)"                    # 通用格式：任何字母和下划线，直到第一个点
    ]
    
    prefix = "unknown"
    for pattern in prefix_patterns
        prefix_match = match(pattern, jobname)
        if prefix_match !== nothing
            prefix = prefix_match.captures[1]
            break
        end
    end
    
    # 定义所有可能的参数匹配模式
    param_patterns = [
        (:b, r"b(\d+\.\d+)"),           # 温度参数b
        (:U, r"U(-?\d+\.\d+)"),         # 相互作用强度U
        (:gw, r"gw(-?\d+\.\d+)"),       # Gutzwiller参数gw
        (:L, r"L(\d+)"),                # 晶格尺寸L
        (:dtau, r"dtau(\d+\.\d+)"),     # 时间步长dtau
        (:lprojgw, r"lprojgw([TF])"),   # 求解投影参数lprojgw
        (:Mz_AFM, r"Mz_AFM(-?\d+\.\d+)"), # AFM参数 Mz_AFM
        (:M_CDW, r"M_CDW(-?\d+\.\d+)"),   # CDW参数M_CDW
        (:xmag, r"xmag(-?\d+\.\d+)")      # 磁场参数xmag
    ]
    
    # 找到所有匹配项并按在原字符串中的位置排序
    all_matches = []
    for (param_key, pattern) in param_patterns
        param_match = match(pattern, jobname)
        if param_match !== nothing
            value_str = param_match.captures[1]
            value, digits = _parse_parameter_value(value_str)
            
            # 跳过NaN值和无效值
            if param_key == :L && value == -1
                continue
            elseif value isa AbstractFloat && isnan(value)
                continue
            elseif param_key == :lprojgw && value === nothing
                continue
            end
            
            # 记录匹配位置以保持原始顺序
            match_pos = param_match.offset
            push!(all_matches, (match_pos, string(param_key), value, digits))
        end
    end
    
    # 按原始位置排序
    sort!(all_matches, by=x->x[1])
    
    # 返回带精度信息的格式，保持原始顺序
    paras_info = [(key, value, digits) for (_, key, value, digits) in all_matches]
    return prefix, paras_info
end

"""
    _parse_parameter_value(value_str::AbstractString) -> (value::Any, digits::Int)

Parse a parameter value string and return the typed value and digits count.
Supports Bool (T/F), Float, and Int values following the design philosophy.
"""
function _parse_parameter_value(value_str::AbstractString)
    value_str = String(value_str)  # Convert to String to handle SubString
    
    # Check for Bool values (T/F format following design philosophy)
    if value_str == "T" || value_str == "true"
        return true, 0
    elseif value_str == "F" || value_str == "false"
        return false, 0
    elseif occursin(".", value_str)
        # Float value
        value = parse(Float64, value_str)
        digits = length(split(value_str, ".")[2])
        return value, digits
    else
        # Integer value
        value = parse(Int64, value_str)
        return value, 0
    end
end

"""
    migrate_legacy_to_new(legacy_jobname::AbstractString) -> String

将legacy格式的作业名转换为新格式，保持参数顺序和精度。

# 参数
- `legacy_jobname::AbstractString`: legacy格式的作业名

# 返回值
- `String`: 新格式的作业名

# 示例
```julia
new_name = migrate_legacy_to_new("proj_fft_honeycomb_exact.b8.000.U4.00.gw0.62.lprojgwF.L9.dtau0.05")
# 返回: "proj_fft_honeycomb_exact_b=8.000_U=4.00_gw=0.62_lprojgw=F_L=9_dtau=0.05"
```
"""
function migrate_legacy_to_new(legacy_jobname::AbstractString)
    # 解析legacy格式
    prefix, params = parse_jobname_legacy(legacy_jobname)
    
    # 生成新格式
    return generate_jobname(prefix, params)
end

"""
    verify_migration(legacy_jobname::AbstractString, new_jobname::AbstractString) -> Bool

验证迁移的正确性，检查参数值和精度是否一致。

# 参数
- `legacy_jobname::AbstractString`: 原始legacy格式作业名
- `new_jobname::AbstractString`: 转换后的新格式作业名

# 返回值
- `Bool`: 如果迁移正确返回true，否则返回false
"""
function verify_migration(legacy_jobname::AbstractString, new_jobname::AbstractString)
    # 解析两种格式
    old_prefix, old_params = parse_jobname_legacy(legacy_jobname)
    new_prefix, new_params = parse_jobname(new_jobname)
    
    # 检查前缀
    if old_prefix != new_prefix
        return false
    end
    
    # 检查参数数量
    if length(old_params) != length(new_params)
        return false
    end
    
    # 检查每个参数
    for ((old_key, old_value, old_digits), (new_key, new_value, new_digits)) in zip(old_params, new_params)
        if old_key != new_key || old_value != new_value || old_digits != new_digits
            return false
        end
    end
    
    return true
end