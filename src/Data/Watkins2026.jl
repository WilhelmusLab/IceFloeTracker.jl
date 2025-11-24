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
import DataFrames: DataFrames, DataFrame, DataFrameRow, nrow, eachrow, subset
import Dates: format
import Images: Gray, rawview, channelview, SegmentedImage

@kwdef struct Case
    loader::GitHubLoader
    metadata::AbstractDict
end

function Case(dfr::DataFrameRow, loader::AbstractLoader)
    metadata = Dict(Symbol(k) => v for (k, v) in pairs(dfr))
    metadata[:case_number] = lpad(dfr.case_number, 3, "0")
    metadata[:region] = dfr.region
    metadata[:date] = format(dfr.start_date, "yyyymmdd")
    metadata[:satellite] = dfr.satellite
    metadata[:pixel_scale] = "250m"
    metadata[:image_side_length] = "100km"
    return Case(; loader, metadata)
end

function metadata(case::Case)::AbstractDict
    return case.metadata
end

function metadata(vcs::Vector{Case})::DataFrame
    return DataFrame(metadata.(vcs))
end

struct Dataset
    loader::AbstractLoader
    metadata::DataFrame
end

function _load_metadata(loader, path)
    df = path |> loader |> load |> DataFrame
    return df
end

function Dataset(;
    loader=GitHubLoader(;
        url="https://github.com/danielmwatkins/ice_floe_validation_dataset/",
        ref="b865acc62f223d6ff14a073a297d682c4c034e5d",
        cache_dir="/tmp/Watkins2026",
    ),
    metadata_path="data/validation_dataset/validation_dataset.csv",
)
    return Dataset(loader, _load_metadata(loader, metadata_path))
end

function metadata(ds::Dataset)::DataFrame
    return ds.metadata
end

function cases(ds::Dataset)::Vector{Case}
    vcs = Case.(eachrow(metadata(ds)), Ref(ds.loader))
    return vcs
end

Base.length(ds::Dataset)::Int = length(cases(ds))
Base.iterate(ds::Dataset) = iterate(cases(ds))
Base.iterate(ds::Dataset, state) = iterate(cases(ds), state)
Base.filter(f::Function, ds::Dataset)::Dataset =
    Dataset(ds.loader, filter(row -> f(row), ds.metadata))
subset(ds::Dataset, args...; kwargs...)::Dataset =
    Dataset(ds.loader, subset(ds.metadata, args...; kwargs...))

function modis_truecolor(case::Case; ext="tiff")
    m = case.metadata
    file = "data/modis/truecolor/$(m[:case_number])-$(m[:region])-$(m[:image_side_length])-$(m[:date]).$(m[:satellite]).truecolor.$(m[:pixel_scale]).$(ext)"
    img = file |> case.loader |> load
    return img
end

function modis_falsecolor(case::Case; ext="tiff")
    m = case.metadata
    file = "data/modis/falsecolor/$(m[:case_number])-$(m[:region])-$(m[:image_side_length])-$(m[:date]).$(m[:satellite]).falsecolor.$(m[:pixel_scale]).$(ext)"
    img = file |> case.loader |> load
    return img
end
function modis_landmask(case::Case; ext="tiff")
    m = case.metadata
    file = "data/modis/landmask/$(m[:case_number])-$(m[:region])-$(m[:image_side_length])-$(m[:date]).$(m[:satellite]).landmask.$(m[:pixel_scale]).$(ext)"
    img = file |> case.loader |> load .|> Gray .|> (x -> x .> 0.5) .|> Gray
    return img
end
function modis_cloudfraction(case::Case; ext="tiff")
    m = case.metadata
    file = "data/modis/cloudfraction/$(m[:case_number])-$(m[:region])-$(m[:image_side_length])-$(m[:date]).$(m[:satellite]).cloudfraction.$(m[:pixel_scale]).$(ext)"
    img = file |> case.loader |> load

    return img
end
function masie_landmask(case::Case; ext="tiff")
    m = case.metadata
    file = "data/masie/landmask/$(m[:case_number])-$(m[:region])-$(m[:image_side_length])-$(m[:date]).masie.landmask.$(m[:pixel_scale]).$(ext)"
    img = file |> case.loader |> load |> (x -> x .> 0.5) .|> Gray
    return img
end
function masie_seaice(case::Case; ext="tiff")
    m = case.metadata
    file = "data/masie/seaice/$(m[:case_number])-$(m[:region])-$(m[:image_side_length])-$(m[:date]).masie.seaice.$(m[:pixel_scale]).$(ext)"
    img = file |> case.loader |> load |> (x -> x .> 0.5) .|> Gray

    return img
end
function validated_binary_floes(case::Case)
    m = case.metadata
    file = "data/validation_dataset/binary_floes/$(m[:case_number])-$(m[:region])-$(m[:date])-$(m[:satellite])-binary_floes.png"
    img = file |> case.loader |> load .|> Gray |> (x -> x .> 0.5) .|> Gray

    return img
end
function validated_labeled_floes(case::Case; ext="tiff")
    m = case.metadata
    file = "data/validation_dataset/labeled_floes/$(m[:case_number])-$(m[:region])-$(m[:date])-$(m[:satellite])-labeled_floes.$(ext)"
    labels = file |> case.loader |> load .|> Int
    img = SegmentedImage(modis_truecolor(case), labels)
    return img
end
function validated_floe_properties(case::Case)::DataFrame
    m = case.metadata
    file = "data/validation_dataset/property_tables/$(m[:satellite])/$(m[:case_number])-$(m[:region])-$(m[:date])-$(m[:satellite])-floe_properties.csv"
    img = file |> case.loader |> load |> DataFrame

    return img
end

end
