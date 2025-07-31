using Downloads: download, RequestError
using CSV
using Dates
using Images: SegmentedImage, Colorant, Gray
using FileIO: load
using CSVFiles
using DataFrames

@kwdef struct ValidationDataCase
    metadata::Union{AbstractDict,Missing} = missing

    # TODO: add types
    modis_truecolor::Union{AbstractArray,Missing} = missing
    modis_falsecolor::Union{AbstractArray,Missing} = missing
    modis_landmask::Union{AbstractArray,Missing} = missing
    modis_cloudfraction::Union{AbstractArray,Missing} = missing
    maisie_landmask::Union{AbstractArray,Missing} = missing
    maisie_seaice::Union{AbstractArray,Missing} = missing

    validated_binary_floes::Union{AbstractArray{Gray{Bool}},Missing} = missing
    validated_labeled_floes::Union{SegmentedImage,Missing} = missing
    validated_floe_properties::Union{CSVFiles.CSVFile,Missing} = missing
end

abstract type ValidationDataLoader end

@kwdef struct Watkins2025GitHub <: ValidationDataLoader
    url::AbstractString = "https://raw.githubusercontent.com/danielmwatkins/ice_floe_validation_dataset/refs/heads/"
    ref::AbstractString = "main"
    dataset_metadata_path::AbstractString = "data/validation_dataset/validation_dataset.csv"
    target_directory::AbstractString = "./Watkins2025GitHub"
end

function (p::ValidationDataLoader)(; kwargs...)
    metadata_url = joinpath(p.url, p.ref, p.dataset_metadata_path)
    metadata_path = joinpath(
        p.target_directory, p.ref, splitpath(p.dataset_metadata_path)[end]
    )
    mkpath(dirname(metadata_path))
    isfile(metadata_path) || download(metadata_url, metadata_path) # Only download if the file doesn't already exist
    metadata = CSV.File(metadata_path)
    @show metadata

    return metadata
end

function load_case(case::CSV.Row, p::Watkins2025GitHub)
    validation_data_dict = Dict()
    validation_data_dict[:metadata] = Dict(
        symbol => case[symbol] for symbol in propertynames(case)
    )

    case_number = lpad(case.case_number, 3, "0")
    region = case.region
    date = Dates.format(case.start_date, "yyyymmdd")
    satellite = case.satellite
    pixel_scale = "250m"
    image_side_length = "100km"
    ext = "tiff"

    output_directory = joinpath(
        p.target_directory,
        p.ref,
        "$(case_number)-$(region)-$(image_side_length)-$(date)-$(satellite)-$(pixel_scale)",
    )
    mkpath(output_directory)

    metadata_path = joinpath(output_directory, "case_metadata.csv")
    isfile(metadata_path) || CSV.write(metadata_path, [case])

    modis_truecolor = (;
        source="data/modis/truecolor/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).truecolor.$(pixel_scale).$(ext)",
        target="modis_truecolor.$(ext)",
        name=:modis_truecolor,
    )
    modis_falsecolor = (;
        source="data/modis/falsecolor/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).falsecolor.$(pixel_scale).$(ext)",
        target="modis_falsecolor.$(ext)",
        name=:modis_falsecolor,
    )
    modis_landmask = (;
        source="data/modis/landmask/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).landmask.$(pixel_scale).$(ext)",
        target="modis_landmask.$(ext)",
        name=:modis_landmask,
    )
    modis_cloudfraction = (;
        source="data/modis/cloudfraction/$(case_number)-$(region)-$(image_side_length)-$(date).$(satellite).cloudfraction.$(pixel_scale).$(ext)",
        target="modis_cloudfraction.$(ext)",
        name=:modis_cloudfraction,
    )
    maisie_landmask = (;
        source="data/masie/landmask/$(case_number)-$(region)-$(image_side_length)-$(date).masie.landmask.$(pixel_scale).$(ext)",
        target="maisie_landmask.$(ext)",
        name=:maisie_landmask,
    )
    maisie_seaice = (;
        source="data/masie/seaice/$(case_number)-$(region)-$(image_side_length)-$(date).masie.seaice.$(pixel_scale).$(ext)",
        target="maisie_seaice.$(ext)",
        name=:maisie_seaice,
    )
    validated_binary_floes = (;
        source="data/validation_dataset/binary_floes/$(case_number)-$(region)-$(date)-$(satellite)-binary_floes.png",
        target="validated_binary_floes.png",
        name=:validated_binary_floes,
    )
    validated_labeled_floes = (;
        source="data/validation_dataset/labeled_floes/$(case_number)-$(region)-$(date)-$(satellite)-labeled_floes.$(ext)",
        target="validated_labeled_floes.$(ext)",
        name=:validated_labeled_floes,
    )
    validated_floe_properties = (;
        source="data/validation_dataset/property_tables/$(satellite)/$(case_number)-$(region)-$(date)-$(satellite)-floe_properties.csv",
        target="validated_floe_properties.csv",
        name=:validated_floe_properties,
    )

    for file_information in [
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
        try
            file_url = joinpath(p.url, p.ref, file_information.source)
            file_path = joinpath(output_directory, file_information.target)
            validation_data_dict[file_information.name] = get_file(file_url, file_path)
        catch e
            if isa(e, RequestError)
                @show "$(file_url) missing"
            else
                rethrow(e)
            end
        end
    end

    # Conversions
    validation_data_dict[:validated_labeled_floes] = SegmentedImage(
        validation_data_dict[:modis_truecolor],
        Int.(validation_data_dict[:validated_labeled_floes]),
    )
    validation_data_dict[:validated_binary_floes] =
        Gray.(Gray.(validation_data_dict[:validated_binary_floes]) .> 0.5)

    validation_data = ValidationDataCase(; validation_data_dict...)
    return validation_data
end

function get_file(file_url, file_path)
    try
        isfile(file_path) || download(file_url, file_path)
        file = load(file_path)
        return file
    catch e
        if isa(e, RequestError)
            @show "$(file_url) missing"
            return nothing
        else
            rethrow(e)
        end
    end
end