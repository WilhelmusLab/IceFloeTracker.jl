module Data

export ValidationDataSet, ValidationDataCase, ValidationDataLoader, Watkins2025GitHub

import Downloads: download, RequestError
import Dates: format
import Images: SegmentedImage, Colorant, Gray
import FileIO: load, save
import CSVFiles: CSVFile
import DataFrames: DataFrame

@kwdef struct ValidationDataSet
    data::Base.Generator
    metadata::DataFrame
end
Base.iterate(iter::ValidationDataSet) = iterate(iter.data)
Base.iterate(iter::ValidationDataSet, state) = iterate(iter.data, state)

@kwdef struct ValidationDataCase
    name::AbstractString = nothing
    metadata::Union{AbstractDict,Nothing} = nothing

    modis_truecolor::Union{AbstractArray,Nothing} = nothing
    modis_falsecolor::Union{AbstractArray,Nothing} = nothing
    modis_landmask::Union{AbstractArray,Nothing} = nothing
    modis_cloudfraction::Union{AbstractArray,Nothing} = nothing
    masie_landmask::Union{AbstractArray,Nothing} = nothing
    masie_seaice::Union{AbstractArray,Nothing} = nothing

    validated_binary_floes::Union{AbstractArray{Gray{Bool}},Nothing} = nothing
    validated_labeled_floes::Union{SegmentedImage,Nothing} = nothing
    validated_floe_properties::Union{CSVFile,Nothing} = nothing
end

"""
Loader for validated ice floe data.
"""
abstract type ValidationDataLoader end

"""

    Watkins2025GitHub(; ref)()
    Watkins2025GitHub(; ref, [url, dataset_metadata_path, cache_dir])(; [case_filter])

Loader for validated ice floe data from [the Watkins 2025 Ice Floe Validation Dataset](https://github.com/danielmwatkins/ice_floe_validation_dataset).

The loader is initialized with a specific `git` ref (tag or commit ID) from which to load the data.

```jldoctest Watkins2025GitHub; setup = :(using IceFloeTracker)
julia> data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
```

`Watkins2025GitHub` fields: 
- `ref`: `git` ref of the commit from which to load the data
- `cache_dir` (optional): local path where the data will be stored, which defaults to `/tmp/Watkins2025/`.
- `url` (optional): URL of the GitHub repository with the dataset
- `dataset_metadata_path` (optional): path within the repository to a CSV file describing the data

The loader is then called with an optional `case_filter` function to filter which cases to load.
This checks each case's metadata, and if the function returns `true`, the case is included in the returned dataset.

```jldoctest Watkins2025GitHub
julia> dataset = data_loader(;case_filter=c -> (
                        c.visible_floes == "yes" &&
                        c.cloud_category_manual == "none" &&
                        c.artifacts == "no"
                    ));
```

The returned `dataset` (a `ValidationDataSet`) has a `metadata` field with a DataFrame of the cases which passed the filter:

```jldoctest Watkins2025GitHub
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
```

The `dataset` contains `ValidationDataCase` objects.
Each ValidationDataCase has metadata fields including:
- `name`: name of the case
- `metadata`: dictionary of metadata for the case, corresponding to a row in the `dataset.metadata` `DataFrame`

Each ValidationDataCase also has data fields including:
- `modis_truecolor`: MODIS true color image
- `modis_falsecolor`: MODIS false color image
- `modis_landmask`: MODIS landmask image
- `modis_cloudfraction`: MODIS cloud fraction image
- `masie_landmask`: MASIE landmask image
- `masie_seaice`: MASIE sea ice image

A ValidationDataCase may have validated data fields including:
- `validated_binary_floes`: binary image of validated floes
- `validated_labeled_floes`: labeled image of validated floes
- `validated_floe_properties`: CSV file of validated floe properties

The `dataset` can be iterated over to get each `ValidationDataCase`:
Example:
```jldoctest Watkins2025GitHub
julia> for case in dataset
           println(case.name * 
                   ": sea ice fraction: " * string(case.metadata[:sea_ice_fraction]) *
                   ", true color image size: " * string(size(case.modis_truecolor)))
       end
011-baffin_bay-100km-20110702-aqua-250m: sea ice fraction: 0.8, true color image size: (400, 400)
014-baffin_bay-100km-20220706-terra-250m: sea ice fraction: 1.0, true color image size: (400, 400)
048-beaufort_sea-100km-20210427-aqua-250m: sea ice fraction: 1.0, true color image size: (400, 400)
048-beaufort_sea-100km-20210427-terra-250m: sea ice fraction: 1.0, true color image size: (400, 400)
054-beaufort_sea-100km-20150516-aqua-250m: sea ice fraction: 1.0, true color image size: (400, 400)
054-beaufort_sea-100km-20150516-terra-250m: sea ice fraction: 1.0, true color image size: (400, 400)
128-hudson_bay-100km-20190415-aqua-250m: sea ice fraction: 1.0, true color image size: (400, 400)
166-laptev_sea-100km-20160904-aqua-250m: sea ice fraction: 1.0, true color image size: (400, 400)
```

!!! info "`dataset` and `dataset.data`" 
    Iterating over the `dataset` is the same as iterating over `dataset.data`, 
    so you could also write `for case in dataset.data...`.)

To get the first case in the dataset, you can use `first(...)`:

```jldoctest Watkins2025GitHub
julia> first(dataset)
ValidationDataCase("011-baffin_bay-100km-20110702-aqua-250m", Dict{Symbol, Any}(:sea_ice_fraction => 0.8, :vi...

julia> first(dataset).validated_labeled_floes
Segmented Image with:
labels map: 400×400 Matrix{Int64}
number of labels: 105

julia> first(dataset).modis_truecolor
400×400 Array{RGBA{N0f8},2} with eltype ColorTypes.RGBA{FixedPointNumbers.N0f8}:
 RGBA{N0f8}(0.094,0.133,0.169,1.0)  RGBA{N0f8}(0.051,0.094,0.118,1.0) ...

```

!!! tip "Cacheing"
    Data are downloaded to the `<cache_dir>/<ref>`, e.g. `/tmp/Watkins2025/a451cd5e62a10309a9640fbbe6b32a236fcebc70/`. 
    If a file of the correct name already exists in that path, if loaded again the cached data will be returned.

    There are no checks to ensure that the cached data are up-to-date, 
    so if the data change in the source for that `ref`, the loader won't load the new data.
    In this case, you can clear the cache by deleting the cache directory, 
    e.g. `rm -r /tmp/Watkins2025/a451cd5e62a10309a9640fbbe6b32a236fcebc70/`.

```
"""
@kwdef struct Watkins2025GitHub <: ValidationDataLoader
    ref::AbstractString
    url::AbstractString = "https://github.com/danielmwatkins/ice_floe_validation_dataset/"
    dataset_metadata_path::AbstractString = "data/validation_dataset/validation_dataset.csv"
    cache_dir::AbstractString = "/tmp/Watkins2025"
end

function (p::ValidationDataLoader)(; case_filter::Function=(case) -> true)
    # Load the metadata
    all_metadata = _load_metadata(p)

    # Filter the metadata
    filtered_metadata = filter(case_filter, all_metadata)

    # Load the data for the filtered metadata
    filtered_data = (_load_case(case, p) for case in eachrow(filtered_metadata))

    return ValidationDataSet(; data=filtered_data, metadata=filtered_metadata)
end

function (p::ValidationDataLoader)(case_filter::Function)
    return p(; case_filter)
end

function _load_metadata(p::ValidationDataLoader)::DataFrame
    metadata_url = joinpath(p.url, "raw", p.ref, p.dataset_metadata_path)
    metadata_path = joinpath(p.cache_dir, p.ref, splitpath(p.dataset_metadata_path)[end])
    mkpath(dirname(metadata_path))

    # Load the metadata file
    isfile(metadata_path) || download(metadata_url, metadata_path) # Only download if the file doesn't already exist
    metadata = DataFrame(load(metadata_path))
    return metadata
end

Base.length(p::ValidationDataLoader) = nrow(_load_metadata(p))

function _load_case(case, p::Watkins2025GitHub)::ValidationDataCase
    data_dict = Dict()
    data_dict[:metadata] = Dict(symbol => case[symbol] for symbol in propertynames(case))

    case_number = lpad(case.case_number, 3, "0")
    region = case.region
    date = format(case.start_date, "yyyymmdd")
    satellite = case.satellite
    pixel_scale = "250m"
    image_side_length = "100km"
    ext = "tiff"

    name = "$(case_number)-$(region)-$(image_side_length)-$(date)-$(satellite)-$(pixel_scale)"
    data_dict[:name] = name

    output_directory = joinpath(p.cache_dir, p.ref, name)
    mkpath(output_directory)

    metadata_path = joinpath(output_directory, "case_metadata.csv")
    isfile(metadata_path) || save(metadata_path, DataFrame(case))

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
    masie_landmask = (;
        source="data/masie/landmask/$(case_number)-$(region)-$(image_side_length)-$(date).masie.landmask.$(pixel_scale).$(ext)",
        target="masie_landmask.$(ext)",
        name=:masie_landmask,
    )
    masie_seaice = (;
        source="data/masie/seaice/$(case_number)-$(region)-$(image_side_length)-$(date).masie.seaice.$(pixel_scale).$(ext)",
        target="masie_seaice.$(ext)",
        name=:masie_seaice,
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
        masie_landmask,
        masie_seaice,
        validated_binary_floes,
        validated_labeled_floes,
        validated_floe_properties,
    ]
        file_url = joinpath(p.url, "raw", p.ref, file_information.source)
        file_path = joinpath(output_directory, file_information.target)
        data_dict[file_information.name] = _get_file(file_url, file_path)
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

function _get_file(file_url, file_path)
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

end
