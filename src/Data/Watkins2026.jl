"""
Module for loading validated ice floe data.
"""
module GitHubData

abstract type AbstractLoader end
import Downloads: RequestError, download

@kwdef struct Loader <: AbstractLoader
    url::AbstractString
    ref::AbstractString
    cache_dir::AbstractString
end

function (p::Loader)(file::AbstractString)::AbstractString
    source = joinpath(p.url, "raw", p.ref, file)
    target = joinpath(p.cache_dir, p.ref, file)
    data = _get_file(source, target)
    return data
end

function _get_file(file_url::AbstractString, file_path::AbstractString)::AbstractString
    @debug "looking for file at $(file_path). File exists: $(isfile(file_path))"
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
    loader::AbstractLoader
    metadata::DataFrameRow
end

function metadata(case::Case)::DataFrameRow
    return case.metadata
end

function loader(case::Case)::AbstractLoader
    return case.loader
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

Base.length(ds::Dataset)::Int = length(cases(ds))
Base.iterate(ds::Dataset) = iterate(cases(ds))
Base.iterate(ds::Dataset, state) = iterate(cases(ds), state)
Base.filter(f::Function, ds::Dataset)::Dataset =
    Dataset(ds.loader, filter(row -> f(row), ds.metadata))
Base.getindex(ds::Dataset, i::Int)::Case = Case(ds.loader, ds.metadata[i, :])
subset(ds::Dataset, args...; kwargs...)::Dataset =
    Dataset(ds.loader, subset(ds.metadata, args...; kwargs...))

function cases(ds::Dataset)::Vector{Case}
    vcs = Case.(Ref(ds.loader), eachrow(metadata(ds)))
    return vcs
end

function metadata(vcs::Vector{Case})::DataFrame
    return DataFrame(metadata.(vcs))
end

function modis_truecolor(case::Case; ext="tiff")
    (; case_number, region, date, satellite, pixel_scale, image_scale) = _filename_parts(
        case
    )
    file = "data/modis/truecolor/$(case_number)-$(region)-$(image_scale)-$(date).$(satellite).truecolor.$(pixel_scale).$(ext)"
    img = file |> case.loader |> load
    return img
end

function modis_falsecolor(case::Case; ext="tiff")
    (; case_number, region, date, satellite, pixel_scale, image_scale) = _filename_parts(
        case
    )
    file = "data/modis/falsecolor/$(case_number)-$(region)-$(image_scale)-$(date).$(satellite).falsecolor.$(pixel_scale).$(ext)"
    img = file |> case.loader |> load
    return img
end

function modis_landmask(case::Case; ext="tiff")
    (; case_number, region, date, satellite, pixel_scale, image_scale) = _filename_parts(
        case
    )
    file = "data/modis/landmask/$(case_number)-$(region)-$(image_scale)-$(date).$(satellite).landmask.$(pixel_scale).$(ext)"
    img = file |> case.loader |> load .|> Gray .|> (x -> x .> 0.1) .|> Gray
    return img
end

function modis_cloudfraction(case::Case; ext="tiff")
    (; case_number, region, date, satellite, pixel_scale, image_scale) = _filename_parts(
        case
    )
    file = "data/modis/cloudfraction/$(case_number)-$(region)-$(image_scale)-$(date).$(satellite).cloudfraction.$(pixel_scale).$(ext)"
    img = file |> case.loader |> load
    return img
end

function validated_binary_floes(case::Case)
    (; case_number, region, date, satellite) = _filename_parts(case)
    file = "data/validation_dataset/binary_floes/$(case_number)-$(region)-$(date)-$(satellite)-binary_floes.png"
    img = file |> case.loader |> load .|> Gray |> (x -> x .> 0.5) .|> Gray
    return img
end

function validated_labeled_floes(case::Case; ext="tiff")
    (; case_number, region, date, satellite) = _filename_parts(case)
    file = "data/validation_dataset/labeled_floes/$(case_number)-$(region)-$(date)-$(satellite)-labeled_floes.$(ext)"
    labels = file |> case.loader |> load .|> Int
    img = SegmentedImage(modis_truecolor(case), labels)
    return img
end

function validated_floe_properties(case::Case)::DataFrame
    (; case_number, region, date, satellite) = _filename_parts(case)
    file = "data/validation_dataset/property_tables/$(satellite)/$(case_number)-$(region)-$(date)-$(satellite)-floe_properties.csv"
    img = file |> case.loader |> load |> DataFrame
    return img
end

function masie_landmask(case::Case; ext="tiff")
    @warn "MASIE landmask data is all zeroes."
    (; case_number, region, date, pixel_scale, image_scale) = _filename_parts(case)
    file = "data/masie/landmask/$(case_number)-$(region)-$(image_scale)-$(date).masie.landmask.$(pixel_scale).$(ext)"
    img = file |> case.loader |> load
    return img
end

function masie_seaice(case::Case; ext="tiff")
    @warn "MASIE sea ice data is all zeroes."
    (; case_number, region, date, pixel_scale, image_scale) = _filename_parts(case)
    file = "data/masie/seaice/$(case_number)-$(region)-$(image_scale)-$(date).masie.seaice.$(pixel_scale).$(ext)"
    img = file |> case.loader |> load
    return img
end

function _filename_parts(case::Case)
    m = metadata(case)
    case_number = lpad(m.case_number, 3, "0")
    region = m.region
    date = format(m.start_date, "yyyymmdd")
    satellite = m.satellite
    pixel_scale = "250m"
    image_scale = "100km"
    return (; case_number, region, date, satellite, pixel_scale, image_scale)
end

end
