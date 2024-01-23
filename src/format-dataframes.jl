using DelimitedFiles
using Printf
using LaTeXStrings
using DataFrames
using TableMetadataTools
using CSV

export format_energy, format_pair_onsite_edge

"""
    format_energy(data_dir::String)
Format the `energy.bin` in `data_dir` into a dataframe and store it in a .csv and .toml file.
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
Format the `pair_onsite_edge.bin`, `pair_onsite_interedges.bin` and `pair_onsite_intraedges.bin` in `data_dir` into a dataframe and store it in a .csv and .toml file.
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