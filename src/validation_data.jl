using Downloads: download, RequestError
using Dates
using Images: SegmentedImage, Colorant, Gray
using FileIO: load
using CSVFiles
using DataFrames

@kwdef struct ValidationDataCase
    name::AbstractString = nothing
    metadata::Union{AbstractDict,Nothing} = nothing

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

"""
Loader for validated ice floe data.
"""
abstract type ValidationDataLoader end

"""
Loader for validated ice floe data structured like https://github.com/danielmwatkins/ice_floe_validation_dataset.

    Struct fields:
    - `url`: URL of the GitHub repository with the dataset
    - `ref`: `git` ref of the commit from which to load the data
    - `dataset_metadata_path`: path within the repository to a CSV file describing the data
    - `cache_dir`: local path where the data will be stored.

    
    Cacheing: 
    Data are downloaded to the `ref` subdirectory of `cache_dir`, e.g. /tmp/Watkins2025/main`. 
    If a file of the correct name already exists in that path, if loaded again the cached data will be returned.
    If the data change in the source for that ref, the loader won't load the new data.
    In that case, it's necessary to delete the cached file.
    A better choice is to use a specific revision `ref`: either a tag, or a commit ID.

    Function arguments:
    - `case_filter`: function run on each metadata entry; 
      if it returns true, then the data from that case is returned 

    Function returns a named tuple with these fields:
    - `metadata`: DataFrame of the cases which passed the `case_filter`
    - `data`: Generator which returns a `ValidationDataCase` for each case which passed the `case_filter`
        
    Example:
    ```
    julia> data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
    julia> dataset = data_loader(;case_filter=c -> (
                            c.visible_floes == "yes" &&
                            c.cloud_category_manual == "none" &&
                            c.artifacts == "no"
                        ));
    julia> dataset.metadata
    8×28 DataFrame
     Row │ case_number  region        start_date  center_lon  center_lat  center_x  center_y  month  sea_ice_fr ⋯
         │ Int64        String31      Date        Float64     Float64     Int64     Int64     Int64  Float64    ⋯
    ─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────
       1 │          11  baffin_bay    2011-07-02    -70.7347     72.3303   -837500  -1737500      7             ⋯
       2 │          14  baffin_bay    2022-07-06    -69.0755     72.3157   -787500  -1762500      7
       3 │          48  beaufort_sea  2021-04-27   -140.612      70.1346  -2162500    212500      4
       4 │          48  beaufort_sea  2021-04-27   -140.612      70.1346  -2162500    212500      4
       5 │          54  beaufort_sea  2015-05-16   -136.675      70.4441  -2137500     62500      5             ⋯
       6 │          54  beaufort_sea  2015-05-16   -136.675      70.4441  -2137500     62500      5
       7 │         128  hudson_bay    2019-04-15    -91.9847     57.853   -2612500  -2437500      4
       8 │         166  laptev_sea    2016-09-04    136.931      79.7507    -37500   1112500      9
                                                                                            20 columns omitted

    julia> first(dataset.data)
    ValidationDataCase("011-baffin_bay-100km-20110702-aqua-250m", Dict{Symbol, Any}(:sea_ice_fraction => 0.8, :vi...

    julia> first(dataset.data).validated_labeled_floes
    Segmented Image with:
    labels map: 400×400 Matrix{Int64}
    number of labels: 105
    ```
"""
@kwdef struct Watkins2025GitHub <: ValidationDataLoader
    ref::AbstractString
    url::AbstractString = "https://github.com/danielmwatkins/ice_floe_validation_dataset/"
    dataset_metadata_path::AbstractString = "data/validation_dataset/validation_dataset.csv"
    cache_dir::AbstractString = "/tmp/Watkins2025"
end

function (p::ValidationDataLoader)(; case_filter::Function=(case) -> true)
    metadata_url = joinpath(p.url, "raw", p.ref, p.dataset_metadata_path)
    metadata_path = joinpath(p.cache_dir, p.ref, splitpath(p.dataset_metadata_path)[end])
    mkpath(dirname(metadata_path))

    # Load the metadata file
    isfile(metadata_path) || download(metadata_url, metadata_path) # Only download if the file doesn't already exist
    all_metadata = DataFrame(load(metadata_path))

    # Filter the metadata
    filtered_metadata = filter(case_filter, all_metadata)

    # Load the data for the filtered metadata
    filtered_data = (_load_case(case, p) for case in eachrow(filtered_metadata))

    return (; data=filtered_data, metadata=filtered_metadata)
end

function _load_case(case, p::Watkins2025GitHub)::ValidationDataCase
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
    isfile(metadata_path) || save(metadata_path, [case])

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
        file_url = joinpath(p.url, "raw", p.ref, file_information.source)
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