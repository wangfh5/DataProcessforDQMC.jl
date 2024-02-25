using DelimitedFiles
using LaTeXStrings
using DataFrames
using TableMetadataTools
using CSV

export format_energy, format_pair_onsite_edge
export format_pair_onsite_bulk, format_pair_onsite_r
export format_onecol, format_rdata, format_kdata

obs_metadata_dict = Dict(
    "pair_onsite_edge" => Dict(
        "label" => L"M_2^{\rm edge}",
        "note" => "s-wave onsite pairing structure factor for edge sites"
    ),
    "pair_onsite_interedges" => Dict(
        "label" => L"M_2^{\rm interedges}",
        "note" => "s-wave onsite pairing structure factor for sites at different edges"
    ),
    "pair_onsite_intraedges" => Dict(
        "label" => L"M_2^{\rm intraedges}",
        "note" => "s-wave onsite pairing structure factor for sites at the same edge"
    ),
    "pair_onsite_bulk" => Dict(
        "label" => L"M_2^{\rm bulk}",
        "note" => "s-wave onsite pairing structure factor for all sites"
    ),
    "pair_onsite_r" => Dict(
        "label" => L"C_{\Delta}(\mathbf{r}_i-\mathbf{r}_j)",
        "note" => "s-wave onsite pairing structure factor (depending on the distance between i and j)"
    ),
    "pair_onsite_k" => Dict(
        "label" => L"C_{\Delta}(\mathbf{k})",
        "note" => "s-wave onsite pairing structure factor in momentum space"
    ),
    "nn_r" => Dict(
        "label" => L"C_{\rm nn}(\mathbf{r}_i-\mathbf{r}_j)",
        "note" => "density-density correlation function (depending on the distance between i and j)"
    ),
    "nn_k" => Dict(
        "label" => L"C_{\rm nn}(\mathbf{k})",
        "note" => "density structure factor in momentum space"
    ),
    "n" => Dict(
        "label" => L"\langle n\rangle",
        "note" => "average occupation number (filling factor) for one sublattice"
    ),
    "spsm_r" => Dict(
        "label" => L"C_{\rm spin}(\mathbf{r}_i-\mathbf{r}_j)=\langle S^+_iS^-_j\rangle",
        "note" => "spin-spin correlation function (depending on the distance between i and j)"
    ),
    "spsm_k" => Dict(
        "label" => L"C_{\rm spin}(\mathbf{k})",
        "note" => "spin susceptibility in momentum space"
    )
)

"""
    write_df(df::DataFrame, data_dir::String, dataname::String)
Write `df` into `data_dir/formatted_data/dataname.csv` and `data_dir/formatted_data/dataname.toml`.
"""
function write_df(df::DataFrame, data_dir::String, dataname::String)
    if !isdir(data_dir * "/formatted_data/")
        mkdir(data_dir * "/formatted_data/")
    end
    CSV.write(data_dir * "/formatted_data/$(dataname).csv", df)
    open(data_dir * "/formatted_data/$(dataname).toml", "w") do io
        print(io, meta2toml(df))
    end
end

"""
    format_energy(data_dir::String)
Format the `energy.bin` in `data_dir` into a dataframe and store it in a `.csv` and `.toml` file.
"""
function format_energy(data_dir::String)
    # read the raw data (energy.bin file)
    energy = readdlm(data_dir * "/energy.bin")
    signs = energy[:,1]
    nbin = size(signs,1)
    # create a dataframe
    df = DataFrame(bin = collect(1:nbin), sig = signs, n = energy[:,3]./signs, Ek = energy[:,5]./signs, Edc = energy[:,7]./signs)
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
    write_df(df, data_dir, "energy")
    return signs
end

## Generic formatters

"""
    format_onecol(data_dir::String, sigs::Array{Float64,1}, obsname::String; ifsave::Bool = false)
A generic formatter for one-column data with each row corresponding to a bin.
Format the `obsname.bin` in `data_dir` into a dataframe and store it in a `.csv` and `.toml` file.
In the `.csv` file, the first column is the bin index, and the second column is the value of the observable divided by `sigs`.
"""
function format_onecol(data_dir::String, sigs::Array{Float64,1}, obsname::String; ifsave::Bool = false)
    # read the raw data (energy.bin file)
    obs = readdlm(data_dir * "/$(obsname).bin")
    nbin = size(sigs,1)
    # create a dataframe
    df = DataFrame(:bin => collect(1:nbin), Symbol(obsname) => obs[:,1]./sigs)
    # add metadata
    caption!(df, "$(obsname).bin")
    metadata!(df, "datadir", data_dir)
    if haskey(obs_metadata_dict, obsname)
        label!(df, Symbol(obsname), obs_metadata_dict[obsname]["label"])
        note!(df, Symbol(obsname), obs_metadata_dict[obsname]["note"])
    else
        @warn "No metadata for $(obsname) is found, potentially do not support this measurement."
    end
    ifsave ? write_df(df, data_dir, obsname) : return df
end
function format_onecol(data_dir::String, sigs::Array{Float64,1}, namelist::Array{String,1}; ifsave::Bool = false)
    [format_onecol(data_dir, sigs, name; ifsave = ifsave) for name in namelist]
end

"""
    format_rdata(data_dir::String, sigs::Array{Float64,1}, obsname::String; ifsave::Bool = false)
A generic formatter for r-dependent data with each row corresponding to a bin and a spatial position.
Format the `obsname.bin` in `data_dir` into a dataframe and store it in a `.csv` and `.toml` file.
For now, only support one sublattice system.
In the `.csv` file, the first column is the bin index, the second and third columns are the x and y coordinates, and the fourth column is the value of the observable divided by `sigs`.
"""
function format_rdata(data_dir::String, sigs::Array{Float64,1}, obsname::String; ifsave::Bool = false)
    # read the raw data (energy.bin file)
    obs = readdlm(data_dir * "/$(obsname).bin")
    nbin = size(sigs,1)
    nsites = Int(size(obs,1)/nbin)
    df = DataFrame(
        :x => obs[:,1], :y => obs[:,2], 
        :bin => repeat(collect(1:nbin),inner = nsites), 
        Symbol(obsname) => obs[:,3]./repeat(sigs, inner = nsites))
    sort!(df,[:x, :y])
    # add metadata
    caption!(df, "$(obsname).bin")
    metadata!(df, "datadir", data_dir)
    label!(df, :x, L"x/a")
    note!(df, :x, "x coordinate of unit cell")
    label!(df, :y, L"y/a")
    note!(df, :y, "y coordinate of unit cell")
    if haskey(obs_metadata_dict, obsname)
        label!(df, Symbol(obsname), obs_metadata_dict[obsname]["label"])
        note!(df, Symbol(obsname), obs_metadata_dict[obsname]["note"])
    else
        @warn "No metadata for $(obsname) is found, potentially do not support this measurement."
    end
    ifsave ? write_df(df, data_dir, obsname) : return df
end
function format_rdata(data_dir::String, sigs::Array{Float64,1}, namelist::Array{String,1}; ifsave::Bool = false)
    [format_rdata(data_dir, sigs, name; ifsave = ifsave) for name in namelist]
end

"""
    format_kdata(data_dir::String, sigs::Array{Float64,1}, obsname::String; ifsave::Bool = false)
"""
function format_kdata(data_dir::String, sigs::Array{Float64,1}, obsname::String; ifsave::Bool = false)
    # read the raw data (energy.bin file)
    obs = readdlm(data_dir * "/$(obsname).bin")
    nbin = size(sigs,1)
    nsites = Int(size(obs,1)/nbin)
    df = DataFrame(
        :kx => obs[:,1], :ky => obs[:,2], 
        :bin => repeat(collect(1:nbin),inner = nsites), 
        Symbol(obsname) => obs[:,3]./repeat(sigs, inner = nsites))
    sort!(df,[:kx, :ky])
    # add metadata
    caption!(df, "$(obsname).bin")
    metadata!(df, "datadir", data_dir)
    label!(df, :kx, L"k_x/(2\pi/a)")
    note!(df, :kx, "x component of lattice momentum")
    label!(df, :ky, L"k_y/(2\pi/a)")
    note!(df, :ky, "y component of lattice momentum")
    if haskey(obs_metadata_dict, obsname)
        label!(df, Symbol(obsname), obs_metadata_dict[obsname]["label"])
        note!(df, Symbol(obsname), obs_metadata_dict[obsname]["note"])
    else
        @warn "No metadata for $(obsname) is found, potentially do not support this measurement."
    end
    ifsave ? write_df(df, data_dir, obsname) : return df
end
function format_kdata(data_dir::String, sigs::Array{Float64,1}, namelist::Array{String,1}; ifsave::Bool = false)
    [format_kdata(data_dir, sigs, name; ifsave = ifsave) for name in namelist]
end

## Specific formatters

"""
    format_pair_onsite_edge(data_dir::String, sigs::Array{Float64,1})
Format the `pair_onsite_edge.bin`, `pair_onsite_interedges.bin` and `pair_onsite_intraedges.bin` in `data_dir` into a dataframe and store it in a `.csv` and `.toml` file.
"""
function format_pair_onsite_edge(data_dir::String, sigs::Array{Float64,1})
    # format each observable
    pair_onsite_edge = format_onecol(data_dir, sigs, "pair_onsite_edge")
    pair_onsite_interedges = format_onecol(data_dir, sigs, "pair_onsite_interedges")
    pair_onsite_intraedges = format_onecol(data_dir, sigs, "pair_onsite_intraedges")
    # combine them into a single dataframe according to the bin index
    df = innerjoin(pair_onsite_edge, pair_onsite_interedges, on = :bin)
    df = innerjoin(df, pair_onsite_intraedges, on = :bin)
    # add table-level metadata
    caption!(df, "pair_onsite_edge.bin, pair_onsite_interedges.bin and pair_onsite_intraedges.bin")
    metadata!(df, "datadir", data_dir)
    # storing metadata persistently
    write_df(df, data_dir, "pair_onsite_edge")
end

"""
    format_pair_onsite_r(data_dir::String, sigs::Array{Float64,1}, L::Int)
Format the `pair_onsite_r.bin` in `data_dir` into a dataframe and store it in a `.csv` and `.toml` file.
For now, only support one sublattice system.
This function support old-fashioned output data, i.e., only output the values without the coordinates, so users need to specify the lattice size `L`.
Just a compatibility function for old-fashioned output, and will be deprecated in the future.
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
    write_df(df, data_dir, "pair_onsite_r")
end