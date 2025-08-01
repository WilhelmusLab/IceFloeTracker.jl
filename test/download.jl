using Downloads: download, RequestError
using CSV
using Dates
using Images: SegmentedImage, Colorant, Gray
using FileIO: load
using CSVFiles
using DataFrames

@kwdef struct ValidationDataCase
    name::AbstractString = nothing
    metadata::Union{AbstractDict,Nothing} = nothing

    # TODO: add types
    modis_truecolor::Union{AbstractArray,Nothing} = nothing
    modis_falsecolor::Union{AbstractArray,Nothing} = nothing
    modis_landmask::Union{AbstractArray,Nothing} = nothing
    modis_cloudfraction::Union{AbstractArray,Nothing} = nothing
    maisie_landmask::Union{AbstractArray,Nothing} = nothing
    maisie_seaice::Union{AbstractArray,Nothing} = nothing

    validated_binary_floes::Union{AbstractArray{Gray{Bool}},Nothing} = nothing
    validated_labeled_floes::Union{SegmentedImage,Nothing} = nothing
    validated_floe_properties::Union{CSVFiles.CSVFile,Nothing} = nothing
end

abstract type ValidationDataLoader end

@kwdef struct Watkins2025GitHub <: ValidationDataLoader
    url::AbstractString = "https://github.com/danielmwatkins/ice_floe_validation_dataset/raw/"
    ref::AbstractString = "main"
    dataset_metadata_path::AbstractString = "data/validation_dataset/validation_dataset.csv"
    cache_dir::AbstractString = "./Watkins2025GitHub"
end

function (p::ValidationDataLoader)(; kwargs...)
    metadata_url = joinpath(p.url, p.ref, p.dataset_metadata_path)
    metadata_path = joinpath(p.cache_dir, p.ref, splitpath(p.dataset_metadata_path)[end])
    mkpath(dirname(metadata_path))
    isfile(metadata_path) || download(metadata_url, metadata_path) # Only download if the file doesn't already exist
    metadata = CSV.File(metadata_path)
    data = (load_case(case, p) for case in metadata)
    return (; data, metadata)
end

function load_case(case::CSV.Row, p::Watkins2025GitHub)
    data_dict = Dict()
    data_dict[:metadata] = Dict(symbol => case[symbol] for symbol in propertynames(case))

    case_number = lpad(case.case_number, 3, "0")
    region = case.region
    date = Dates.format(case.start_date, "yyyymmdd")
    satellite = case.satellite
    pixel_scale = "250m"
    image_side_length = "100km"
    ext = "tiff"

    name = "$(case_number)-$(region)-$(image_side_length)-$(date)-$(satellite)-$(pixel_scale)"
    data_dict[:name] = name

    output_directory = joinpath(p.cache_dir, p.ref, name)
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
        file_url = joinpath(p.url, p.ref, file_information.source)
        file_path = joinpath(output_directory, file_information.target)
        data_dict[file_information.name] = get_file(file_url, file_path)
    end

    # Conversions
    if !isnothing(data_dict[:validated_labeled_floes]) &&
        !isnothing(data_dict[:modis_truecolor])
        data_dict[:validated_labeled_floes] = SegmentedImage(
            data_dict[:modis_truecolor], Int.(data_dict[:validated_labeled_floes])
        )
    end
    if !isnothing(data_dict[:validated_binary_floes])
        data_dict[:validated_binary_floes] =
            Gray.(Gray.(data_dict[:validated_binary_floes]) .> 0.5)
    end

    data_struct = ValidationDataCase(; data_dict...)
    return data_struct
end

function get_file(file_url, file_path)
    @debug "looking for file at $(file_path). File exists: $(isfile(file_path))"
    if !isfile(file_path)
        try
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
    file = load(file_path)
    return file
end