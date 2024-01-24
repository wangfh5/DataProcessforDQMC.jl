using DelimitedFiles
using LaTeXStrings
using DataFrames
using TableMetadataTools
using CSV

export format_energy, format_pair_onsite_edge
export format_pair_onsite_bulk, format_pair_onsite_r, format_pair_onsite_k

"""
    format_energy(data_dir::String)
Format the `energy.bin` in `data_dir` into a dataframe and store it in a `.csv` and `.toml` file.
"""
function format_energy(data_dir::String)
    # read the raw data (energy.bin file)
    energy = readdlm(data_dir * "/energy.bin")
    signs = energy[:,1]
    # create a dataframe
    df = DataFrame(sig = signs, n = energy[:,3]./signs, Ek = energy[:,5]./signs, Edc = energy[:,7]./signs)
    # add metadata
    caption!(df, "energy.bin")
    metadata!(df, "datadir", data_dir)
    label!(df, :sig, L"\langle \mathrm{sign} \rangle")
    note!(df, :sig, "average sign")
    label!(df, :n, L"\langle n\rangle")
    note!(df, :n, "average occupation number (filling factor)")
    label!(df, :Ek, L"\langle H_t/t \rangle")
    note!(df, :Ek, "average hopping energy")
    label!(df, :Edc, L"\langle H_{U}/U \rangle")
    note!(df, :Edc, "average double-occupation energy")
    # storing metadata persistently
    if !isdir(data_dir * "/formatted_data/")
        mkdir(data_dir * "/formatted_data/")
    end
    CSV.write(data_dir * "/formatted_data/energy.csv", df)
    open(data_dir * "/formatted_data/energy.toml", "w") do io
        print(io, meta2toml(df))
    end
    return signs
end

"""
    format_pair_onsite_edge(data_dir::String, sigs::Array{Float64,1})
Format the `pair_onsite_edge.bin`, `pair_onsite_interedges.bin` and `pair_onsite_intraedges.bin` in `data_dir` into a dataframe and store it in a `.csv` and `.toml` file.
"""
function format_pair_onsite_edge(data_dir::String, sigs::Array{Float64,1})
    # read the raw data (energy.bin file)
    pair_onsite_edge = readdlm(data_dir * "/pair_onsite_edge.bin")
    pair_onsite_interedges = readdlm(data_dir * "/pair_onsite_interedges.bin")
    pair_onsite_intraedges = readdlm(data_dir * "/pair_onsite_intraedges.bin")
    # create a dataframe
    df = DataFrame(pair_onsite_edge = pair_onsite_edge[:,1]./sigs, pair_onsite_interedges = pair_onsite_interedges[:,1]./sigs, pair_onsite_intraedges = pair_onsite_intraedges[:,1]./sigs)
    # add metadata
    caption!(df, "pair_onsite_edge.bin, pair_onsite_interedges.bin and pair_onsite_intraedges.bin")
    metadata!(df, "datadir", data_dir)
    label!(df, :pair_onsite_edge, L"M_2^{\rm edge}")
    note!(df, :pair_onsite_edge, "s-wave onsite pairing structure factor for edge sites")
    label!(df, :pair_onsite_interedges, L"M_2^{\rm interedges}")
    note!(df, :pair_onsite_interedges, "s-wave onsite pairing structure factor for sites at different edges")
    label!(df, :pair_onsite_intraedges, L"M_2^{\rm intraedges}")
    note!(df, :pair_onsite_intraedges, "s-wave onsite pairing structure factor for sites at the same edge")
    # storing metadata persistently
    if !isdir(data_dir * "/formatted_data/")
        mkdir(data_dir * "/formatted_data/")
    end
    CSV.write(data_dir * "/formatted_data/pair_onsite_edge.csv", df)
    open(data_dir * "/formatted_data/pair_onsite_edge.toml", "w") do io
        print(io, meta2toml(df))
    end
end

"""
    format_pair_onsite_bulk(data_dir::String, sigs::Array{Float64,1})
Format the `pair_onsite_bulk.bin` in `data_dir` into a dataframe and store it in a `.csv` and `.toml` file.
"""
function format_pair_onsite_bulk(data_dir::String, sigs::Array{Float64,1})
    # read the raw data (energy.bin file)
    pair_onsite_bulk = readdlm(data_dir * "/pair_onsite_bulk.bin")
    # create a dataframe
    df = DataFrame(pair_onsite_bulk = pair_onsite_bulk[:,1]./sigs)
    # add metadata
    caption!(df, "pair_onsite_bulk.bin")
    metadata!(df, "datadir", data_dir)
    label!(df, :pair_onsite_bulk, L"M_2^{\rm bulk}")
    note!(df, :pair_onsite_bulk, "s-wave onsite pairing structure factor for all sites")
    # storing metadata persistently
    if !isdir(data_dir * "/formatted_data/")
        mkdir(data_dir * "/formatted_data/")
    end
    CSV.write(data_dir * "/formatted_data/pair_onsite_bulk.csv", df)
    open(data_dir * "/formatted_data/pair_onsite_bulk.toml", "w") do io
        print(io, meta2toml(df))
    end
end

"""
    format_pair_onsite_r(data_dir::String, sigs::Array{Float64,1}, L::Int)
Format the `pair_onsite_r.bin` in `data_dir` into a dataframe and store it in a `.csv` and `.toml` file.
For now, only support one sublattice system.
This function support old-fashioned output data, i.e., only output the values without the coordinates, so users need to specify the lattice size `L`.
"""
function format_pair_onsite_r(data_dir::String, sigs::Array{Float64,1}, L::Int)
    # read the raw data (energy.bin file)
    pair_onsite_r = readdlm(data_dir * "/pair_onsite_r.bin")
    nbin = Int(size(pair_onsite_r, 1)/(L^2))
    # create a dataframe
    df_tmp = DataFrame()
    for i in 1:nbin
        insertcols!(df_tmp, i, Symbol(string(i)) => pair_onsite_r[(i-1)*L^2+1:i*L^2,1]./sigs[i])
    end
    # transpose the table
    insertcols!(df_tmp,1, :imj => collect(1:L^2) )
    df = stack(df_tmp,Not(:imj))
    rename!(df,:variable=>:bin,:value=>:pair_onsite_r)
    df[!, :bin] = parse.(Int, df[!, :bin])
    sort!(df,:imj)
    # add metadata
    caption!(df, "pair_onsite_r.bin")
    metadata!(df, "datadir", data_dir)
    label!(df, :imj, L"\mathbf{r}_i-\mathbf{r}_j")
    note!(df, :imj, "distance between i and j")
    label!(df, :pair_onsite_r, L"C_{\Delta}(\mathbf{r}_i-\mathbf{r}_j)")
    note!(df, :pair_onsite_r, "s-wave onsite pairing structure factor (depending on the distance between i and j)")
    # storing metadata persistently
    if !isdir(data_dir * "/formatted_data/")
        mkdir(data_dir * "/formatted_data/")
    end
    CSV.write(data_dir * "/formatted_data/pair_onsite_r.csv", df)
    open(data_dir * "/formatted_data/pair_onsite_r.toml", "w") do io
        print(io, meta2toml(df))
    end
end

"""
    format_pair_onsite_r(data_dir::String, sigs::Array{Float64,1})
Format the `pair_onsite_r.bin` in `data_dir` into a dataframe and store it in a `.csv` and `.toml` file.
For now, only support one sublattice system.
This function support new-fashioned output data, i.e., output both the coordinates and the values in each line.
"""
function format_pair_onsite_r(data_dir::String, sigs::Array{Float64,1})
    # read the raw data (energy.bin file)
    pair_onsite_r = readdlm(data_dir * "/pair_onsite_r.bin")
    # create a dataframe
    nbin = size(sigs,1)
    nsites = size(pair_onsite_k,1)/nbin
    df = DataFrame(
        imj_x = pair_onsite_r[:,1], imj_y = pair_onsite_r[:,2], 
        bin = repeat(collect(1:nbin),inner = nsites), 
        pair_onsite_r = pair_onsite_r[:,3]./repeat(sigs, inner = nsites))
    sort!(df,[:imj_x, :imj_y])
    # add metadata
    caption!(df, "pair_onsite_r.bin")
    metadata!(df, "datadir", data_dir)
    label!(df, :imj_x, L"x_i-x_j")
    note!(df, :imj_x, "x coordinate difference between i and j")
    label!(df, :imj_y, L"y_i-y_j")
    note!(df, :imj_y, "y coordinate difference between i and j")
    label!(df, :pair_onsite_r, L"C_{\Delta}(\mathbf{r}_i-\mathbf{r}_j)")
    note!(df, :pair_onsite_r, "s-wave onsite pairing structure factor (depending on the distance between i and j)")
    # storing metadata persistently
    if !isdir(data_dir * "/formatted_data/")
        mkdir(data_dir * "/formatted_data/")
    end
    CSV.write(data_dir * "/formatted_data/pair_onsite_r.csv", df)
    open(data_dir * "/formatted_data/pair_onsite_r.toml", "w") do io
        print(io, meta2toml(df))
    end
end

"""
    format_pair_onsite_k(data_dir::String, sigs::Array{Float64,1})
Format the `pair_onsite_k.bin` in `data_dir` into a dataframe and store it in a ``.csv` and `.toml` file.
For now, only support one sublattice system.
"""
function format_pair_onsite_k(data_dir::String, sigs::Array{Float64,1})
    # read the raw data (energy.bin file)
    pair_onsite_k = readdlm(data_dir * "/pair_onsite_k.bin")
    # create a dataframe
    nbin = size(sigs,1)
    nsites = Int(size(pair_onsite_k,1)/nbin)
    df = DataFrame(
        kx = pair_onsite_k[:,1], ky = pair_onsite_k[:,2], 
        bin = repeat(collect(1:nbin),inner = nsites),
        pair_onsite_k = pair_onsite_k[:,3]./repeat(sigs, inner = nsites))
    sort!(df,[:kx, :ky])
    # add metadata
    caption!(df, "pair_onsite_k.bin")
    metadata!(df, "datadir", data_dir)
    label!(df, :kx, L"k_x")
    note!(df, :kx, "x component of momentum")
    label!(df, :ky, L"k_y")
    note!(df, :ky, "y component of momentum")
    label!(df, :pair_onsite_k, L"C_{\Delta}(\mathbf{k})")
    note!(df, :pair_onsite_k, "s-wave onsite pairing structure factor in momentum space")
    # storing metadata persistently
    if !isdir(data_dir * "/formatted_data/")
        mkdir(data_dir * "/formatted_data/")
    end
    CSV.write(data_dir * "/formatted_data/pair_onsite_k.csv", df)
    open(data_dir * "/formatted_data/pair_onsite_k.toml", "w") do io
        print(io, meta2toml(df))
    end
end