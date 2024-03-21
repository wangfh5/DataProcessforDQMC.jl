export read_format_data

para_metadata_dict = Dict(
    "U" => Dict("label" => L"U/t", "note" => "Hubbard interaction strength"),
    "Delta_ppm" => Dict("label" => L"\Delta/t", "note" => "p pm ip pairing strength"),
    "mu" => Dict("label" => L"\mu/t", "note" => "chemical potential"),
    "L" => Dict("label" => L"L", "note" => "system linear size"),
    "beta" => Dict("label" => L"\beta t", "note" => "inverse temperature")
)

function extract_parameters(str::String)
    # str = "proj_majo_cy.b16.000.U-4.00.Delta_ppm0.4.mu-0.5.L8.dtau0.05"
    # 使用正则表达式来提取键值对
    matches = collect(eachmatch(r"([a-zA-Z_]+)-?([\d.-]+)\.?", str))
    # 将提取的键值对存储到字典中
    dict = Dict{String, Number}()
    for match in matches
        key = match.captures[1]
        value = match.captures[2]
        # finetune the key name
        if key == "b"
            key = "beta"
        end
        # 去掉参数值末尾的点
        if endswith(value, ".")
            value = chop(value, tail=1)
        end
        # 不保存参数值为空的键值对
        if value != ""
            if occursin(".", value)
                dict[key] = parse(Float64, value)
            else
                dict[key] = parse(Int, value)
            end
        end
    end
    dict
end

"""
    insert_parameters_cols(df::DataFrame,data_dir::String)
Exact parameters from `data_dir` and insert them as columns into `df`.
e.g. `data_dir = ".../proj_majo_cy.b16.000.U-4.00.Delta_ppm0.4.mu-0.5.L8.dtau0.05/"` 
will extract `beta`, `U`, `Delta_ppm`, `mu`, `L` and `dtau` from it.
"""
@inline function insert_parameters_cols(df::DataFrame,data_dir::String)
    # extract parameters from data_dir
    maindir = splitpath(data_dir)[end]
    paras = extract_parameters(maindir)
    # add parameters column
    # if there are many column, can use array and for loop to add columns; or try join tables
    for (i, (key, value)) in enumerate(paras)
        insertcols!(df, i, Symbol(key) => value)
        if haskey(para_metadata_dict, key)
            label!(df, Symbol(key), para_metadata_dict[key]["label"])
            note!(df, Symbol(key), para_metadata_dict[key]["note"])
        end
    end
    metadata!(df, "npara", length(paras))
    setmetadatastyle!(==("npara"), df)
    df
end

"""
    read_format_data(dataname::String, data_dir::String)
Read `dataname.csv` and `dataname.toml` in `data_dir`.
"""
function read_format_data(dataname::String,data_dir::String;lpara=true)
    # data_dir = "/home/wangfh5/Projects/archieve/majo-ppmip-data/proj_majo_cy.b10.000.U-5.20.Delta_ppm0.4.mu-0.5.L8.dtau0.05/"
    # dataname = "energy"
    # dataname = "pair_onsite_edge"
    tmp_df = CSV.read(data_dir * "formatted_data/$(dataname).csv", DataFrame)
    open(data_dir * "formatted_data/$(dataname).toml", "r") do io
        toml2meta!(tmp_df, io)
    end
    metadata!(tmp_df, "ndata", size(tmp_df, 2)-1)
    setmetadatastyle!(==("ndata"), tmp_df)
    if lpara
        tmp_df = insert_parameters_cols(tmp_df,data_dir)
    end
    tmp_df
end