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

    modis_truecolor = (;
        source="data/modis/truecolor/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).truecolor.$(pixel_scale).$(ext)",
        target="modis_truecolor.$(ext)",
    )
    modis_falsecolor = (;
        source="data/modis/falsecolor/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).falsecolor.$(pixel_scale).$(ext)",
        target="modis_falsecolor.$(ext)",
    )
    modis_landmask = (;
        source="data/modis/landmask/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).landmask.$(pixel_scale).$(ext)",
        target="modis_landmask.$(ext)",
    )
    modis_cloudfraction = (;
        source="data/modis/cloudfraction/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).cloudfraction.$(pixel_scale).$(ext)",
        target="modis_cloudfraction.$(ext)",
    )
    maisie_landmask = (;
        source="data/masie/landmask/$(case_number)-$(region)-$(image_side_length)-$(date).masie.landmask.$(pixel_scale).$(ext)",
        target="maisie_landmask.$(ext)",
    )
    maisie_seaice = (;
        source="data/masie/seaice/$(case_number)-$(region)-$(image_side_length)-$(date).masie.seaice.$(pixel_scale).$(ext)",
        target="maisie_seaice.$(ext)",
    )
    validated_binary_floes = (;
        source="data/validation_dataset/binary_floes/$(case_number)-$(region)-$(date)-$(satellite)-binary_floes.png",
        target="validated_binary_floes.png",
    )
    validated_labeled_floes = (;
        source="data/validation_dataset/labeled_floes/$(case_number)-$(region)-$(date)-$(satellite)-labeled_floes.$(ext)",
        target="validated_labeled_floes.$(ext)",
    )
    validated_floe_properties = (;
        source="data/validation_dataset/property_tables/$(satellite)/$(case_number)-$(region)-$(date)-$(satellite)-floe_properties.csv",
        target="validated_floe_properties.csv",
    )

    output_directory = joinpath(
        p.target_directory,
        "$(case_number)-$(region)-$(image_side_length)-$(date)-$(satellite)-$(pixel_scale)",
    )

    for image in [
        modis_truecolor,
        modis_falsecolor,
        modis_landmask,
        modis_cloudfraction,
        maisie_landmask,
        maisie_seaice,
        validated_binary_floes,
        validated_labeled_floes,
        validated_floe_properties,
    ]
        image_url = joinpath(p.url, image.source)
        image_path = joinpath(output_directory, image.target)
        isfile(image_path) && continue  # don't download a second time if we already have the file
        mkpath(dirname(image_path))
        @debug "downloading $(image_url) to $(image_path)"
        download(image_url, image_path)
    end

    return nothing
end