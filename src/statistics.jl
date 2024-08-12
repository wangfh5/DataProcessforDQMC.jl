using Statistics

"""
    bin_filter(df::DataFrame)
Remove the first `ntherm` bins from `df`. `df` should have a column named `bin`.
"""
@inline function bin_filter(df::DataFrame;ntherm=1::Int64)
    # filter data
    df_filter = df[df.bin .> ntherm,:]
    # remove the bin column
    select!(df_filter, Not(:bin))
    df_filter
end

"""
    k_filter(df::DataFrame, k::Tuple{Float64,Float64})
Filter the data at k = (kx,ky).
"""
@inline function k_filter(df::DataFrame, k::Tuple{Float64,Float64})
    # filter data
    df_filter = df[(df.kx .== k[1]) .& (df.ky .== k[2]),:]
    # remove the kx, ky column
    select!(df_filter, Not(:kx))
    select!(df_filter, Not(:ky))
    df_filter
end

"""
    error(data)
Calculate the error of the data (3σ).
"""
function error(data)
    3*std(data)./sqrt(length(data))
end

"""
    mean_ns(data)
mean without spikes: Calculate the mean of the data after removing the largest and smallest value.
"""
function mean_ns(data)
    sort_data = sort(data)
    mean(sort_data[2:end-1])
end

"""
    error_ns(data)
error without spikes: Calculate the error of the data after removing the largest and smallest value.
"""
function error_ns(data)
    n = length(data)
    sort_data = sort(data)
    3*std(sort_data[2:end-1])/sqrt(n-2)
end

"""
    statistics_columns(df_tmp::DataFrame)
average and error of each column of a dataframe with only `bin` and data columns.
"""
function statistics_columns(df_tmp::DataFrame)
    df_tmp = bin_filter(df_tmp)
    # calculate the mean and error of each column
    mean_values = mean.(eachcol(df_tmp))
    error_values = error.(eachcol(df_tmp))
    # 创建一个存储统计信息的空 DataFrame
    df_stats = DataFrame()
    # 将每一列的平均值和标准差添加到新的 DataFrame 中
    for (colname, mean_val, error_val) in zip(names(df_tmp), mean_values, error_values)
        df_stats[!, Symbol(colname, "_mean")] = [mean_val]
        df_stats[!, Symbol(colname, "_error")] = [error_val]
    end
    return df_stats
end


"""
    statistics_columns_withparas(df_tmp::DataFrame, paras)
average and error of a dataframe with `bin` column, parameter columns specified by `paras` (symbols) and data columns.
"""
function statistics_columns_withparas(df_tmp; paras=nothing)
    df_tmp = bin_filter(df_tmp)
    # find the data names
    datavec = setdiff(Symbol.(names(df_tmp)), paras)
    datamat = reshape(datavec, 1, :)
    # group and combine
    gdf = groupby(df_tmp, paras)
    df_stat = combine(gdf, datamat .=> [mean, error])
    return df_stat
end
