## These function accept the filename of the data file and then do analysis
export RenyiNegativity, RenyiNegativity_all

function RenyiNegativity(filename::String, filedir::String=pwd();printLA=nothing)
    # 打开文件
    filepath = joinpath(filedir, filename)

    # 使用 readdlm 读取整个文件
    # 这里假设文件中的数据以空格分隔
    data = readdlm(filepath, Float64)
    # remove the first two bin
    data = data[3:end,:]
    # remove the max and min
    data = filter(data, 1)

    # deduce L and rank from the data
    L = size(data,2) ÷ 2 - 1
    # example filename: expRenyiN3.bin, expRenyiN4_TW.bin
    rank = parse(Int, split(filename, "expRenyiN")[2][1])
    quantity_name = split(filename, ".bin")[1]

    # 输出结果
    expRenyiN = zeros(L+1,2)
    RenyiN = zeros(L+1,2)
    for i in 1:L+1
        vectmp = data[:,2i-1] # take real part
        if mean(vectmp) > 0.0
            mn = mean(vectmp)
            err = std(vectmp)/sqrt(length(vectmp))
            expRenyiN[i,:] = [mn, err]
            RenyiN[i,:] = [log(mn)/(1-rank), abs(err/mn/(1-rank))]
            if (printLA == nothing) || ((i-1) ∈ printLA)
                println("$(quantity_name) $(i-1) $(expRenyiN[i,:]) $(RenyiN[i,:])")
            end
        end
    end
    # @show expRenyiN
    # @show RenyiN
    return expRenyiN, RenyiN
end

function RenyiNegativity_all(filedir::String=pwd();maxrank::Int=4, kwargs...)
    # find all the files with name expRenyiN*.bin or expRenyiN*_TW.bin, where * is an integer
    filenames = Base.filter(x->occursin(r"expRenyiN\d+\.bin", x) || occursin(r"expRenyiN\d+_TW\.bin", x), readdir(filedir))
    # sort the filenames
    filenames = sort(filenames)
    # get the number of files
    nfiles = length(filenames)
    # for each file, create an array to store the average and error of Renyi negativity
    # and then append the result to an opened JLD2 file
    file = jldopen("RenyiNall.jld2", "w") do file
        for i in 1:nfiles
            filename = filenames[i]
            suffix = split(filename, "expRenyiN")[2][1:end-4]
            rank = parse(Int, suffix[1])
            if rank <= maxrank
                expRenyiN, RenyiN = RenyiNegativity(filename, filedir;kwargs...)
                # save the result to a JLD2 file
                file["expRenyiN$suffix"] = expRenyiN
                file["RenyiN$suffix"] = RenyiN
            end
        end
    end
end

"""
Remove the max and min of the data
"""
function filter(onecolumndata::Array{Float64,1}, dropnum::Int64=1)
    onecolumndata_sort = sort(onecolumndata)
    return onecolumndata_sort[1+dropnum:end-dropnum]
end
function filter(multicolumnsdata::Array{Float64,2}, dropnum::Int64=1)
    return hcat([filter(multicolumnsdata[:,i], dropnum) for i in 1:size(multicolumnsdata,2)]...)
end