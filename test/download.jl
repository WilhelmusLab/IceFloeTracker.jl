using Downloads: download
using CSV
using Dates

abstract type ValidationData end

@kwdef struct Watkins2025GitHub <: ValidationData
    url::AbstractString = "https://raw.githubusercontent.com/danielmwatkins/ice_floe_validation_dataset/refs/heads/main/"
    dataset_metadata_path::AbstractString = "data/validation_dataset/validation_dataset.csv"
    target_directory::AbstractString = "./Watkins2025GitHub"
end

function (p::ValidationData)(; kwargs...)
    mkpath(p.target_directory)

    # TODO: cache this and use the commit SHA to invalidate caches
    metadata_url = joinpath(p.url, p.dataset_metadata_path)
    metadata_path = joinpath(p.target_directory, splitpath(p.dataset_metadata_path)[end])
    download(metadata_url, metadata_path)
    metadata = CSV.File(metadata_path)
    @show metadata

    return metadata
end

function load_case(case::CSV.Row, p::Watkins2025GitHub)
    case_number = lpad(case.case_number, 3, "0")
    region = case.region
    date = Dates.format(case.start_date, "yyyymmdd")
    satellite = case.satellite
    pixel_scale = "250m"
    image_side_length = "100km"
    ext = "tiff"
    modis_truecolor = "data/modis/truecolor/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).truecolor.$(pixel_scale).$(ext)"
    modis_falsecolor = "data/modis/falsecolor/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).falsecolor.$(pixel_scale).$(ext)"
    modis_landmask = "data/modis/landmask/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).landmask.$(pixel_scale).$(ext)"
    modis_cloudfraction = "data/modis/cloudfraction/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).cloudfraction.$(pixel_scale).$(ext)"
    maisie_landmask = "data/masie/landmask/$(case_number)-$(region)-$(image_side_length)-$(date).masie.landmask.$(pixel_scale).$(ext)"
    maisie_seaice = "data/masie/seaice/$(case_number)-$(region)-$(image_side_length)-$(date).masie.seaice.$(pixel_scale).$(ext)"

    for image_path in [
        modis_truecolor,
        modis_falsecolor,
        modis_landmask,
        modis_cloudfraction,
        maisie_landmask,
        maisie_seaice,
    ]
        image_url = joinpath(p.url, image_path)
        image_path = joinpath(p.target_directory, image_path)
        mkpath(dirname(image_path))
        download(image_url, image_path)
    end

    return nothing
end