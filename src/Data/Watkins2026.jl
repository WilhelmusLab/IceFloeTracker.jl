"""
Module for loading validated ice floe data.
"""
module GitHubData

abstract type AbstractLoader end

@kwdef struct Loader <: AbstractLoader
    url::AbstractString
    ref::AbstractString
    cache_dir::AbstractString
end

function (p::Loader)(file::AbstractString)::AbstractString
    @info file
    source = joinpath(p.url, "raw", p.ref, file)
    target = joinpath(p.cache_dir, p.ref, file)
    data = _get_file(source, target)
    return data
end

function _get_file(file_url::AbstractString, file_path::AbstractString)::AbstractString
    @info "looking for file at $(file_path). File exists: $(isfile(file_path))"
    if !isfile(file_path)
        try
            mkpath(dirname(file_path))
            download(file_url, file_path)
        catch e
            if isa(e, RequestError)
                @debug "nothing at $(file_url)"
                return nothing
            else
                rethrow(e)
            end
        end
    end
    return file_path
end

end

module Watkins2026
export Dataset, Case, metadata

import ..GitHubData: Loader as GitHubLoader, AbstractLoader
import FileIO: load, save
import DataFrames: DataFrame, DataFrameRow, nrow, eachrow
import Dates: format

@kwdef struct Case
    loader::GitHubLoader
    metadata::AbstractDict

    case_number::AbstractString
    region::AbstractString
    date::AbstractString
    satellite::AbstractString
    pixel_scale::AbstractString
    image_side_length::AbstractString
end

function Case(dfr::DataFrameRow, loader::AbstractLoader)
    return Case(;
        loader,
        metadata=Dict(Symbol(k) => v for (k, v) in pairs(dfr)),
        case_number=lpad(dfr.case_number, 3, "0"),
        region=dfr.region,
        date=format(dfr.start_date, "yyyymmdd"),
        satellite=dfr.satellite,
        pixel_scale="250m",
        image_side_length="100km",
    )
end

function metadata(case::Case)::AbstractDict
    return case.metadata
end

@kwdef struct Dataset
    loader::AbstractLoader = GitHubLoader(;
        url="https://github.com/danielmwatkins/ice_floe_validation_dataset/",
        ref="b865acc62f223d6ff14a073a297d682c4c034e5d",
        cache_dir="/tmp/Watkins2026",
    )
    dataset_metadata_path::AbstractString = "data/validation_dataset/validation_dataset.csv"
end

function metadata(ds::Dataset)::DataFrame
    csv_path = ds.loader(ds.dataset_metadata_path)
    csv = load(csv_path)
    df = DataFrame(csv)
    return df
end

function cases(ds::Dataset)::Vector{Case}
    df = metadata(ds)
    vcs = Case.(eachrow(df), Ref(ds.loader))
    return vcs
end

Base.length(ds::Dataset)::Int = length(cases(ds))
Base.iterate(ds::Dataset)::Vector{Case} = iterate(cases(ds))
Base.iterate(ds::Dataset, state)::Vector{Case} = iterate(cases(ds), state)

function metadata(vcs::Vector{Case})::DataFrame
    return DataFrame(metadata.(vcs))
end

end
