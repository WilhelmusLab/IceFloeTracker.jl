export Watkins2026Dataset,
    metadata,
    name,
    modis_truecolor,
    modis_falsecolor,
    modis_landmask,
    modis_cloudfraction,
    masie_landmask,
    masie_seaice,
    validated_binary_floes,
    validated_labeled_floes,
    validated_floe_properties

import FileIO: load
using CSVFiles
import DataFrames: DataFrame
import Dates: format
import Images: Gray, SegmentedImage

"""
    Watkins2026Dataset()
    Watkins2026Dataset(; [ref, url, dataset_metadata_path, cache_dir])

Validated ice floe data from [the Watkins 2026 Ice Floe Validation Dataset](https://github.com/danielmwatkins/ice_floe_validation_dataset).

The dataset is initialized with a specific `git` tag, branch or commit ID from which to load the data.

```jldoctest Watkins2026Dataset; setup = :(using IceFloeTracker)
julia> dataset = Watkins2026Dataset(; ref="b865acc62f223d6ff14a073a297d682c4c034e5d")
```

`Watkins2026Dataset` fields: 
- `ref` (optional): `git` tag, commit-id or branch from which to load the data
- `cache_dir` (optional): local path where the data will be stored, which defaults to `/tmp/Watkins2026/`.
- `url` (optional): URL of the GitHub repository with the dataset
- `dataset_metadata_path` (optional): path within the repository to a CSV file describing the data

The dataset can be filtered using `Base.filter` or `DataFrames.subset`.
This checks each case's metadata, and if the function returns `true`, the case is included in the returned dataset.

```jldoctest Watkins2026Dataset
julia> dataset = filter(c -> (
                        c.visible_floes == "yes" &&
                        c.cloud_category_manual == "none" &&
                        c.artifacts == "no"
                    ), dataset);
```

Equivalently:
```jldoctest Watkins2026Dataset
julia> dataset = subset(dataset, 
                        :visible_floes => c -> c .== "yes",
                        :cloud_category_manual => c -> c .== "none",
                        :artifacts => c -> c .== "no",
                    );
```

The returned `dataset` (a `Dataset`) has a `metadata` accessor which returns a DataFrame of the cases which passed the filter:

```jldoctest Watkins2026Dataset
julia> metadata(dataset)
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

The `dataset` contains `Case` objects.
Each Case has metadata fields including:
- `name`: name of the case
- `metadata`: dictionary of metadata for the case, corresponding to a row in the `dataset.metadata` `DataFrame`

Each Case has functions to access its contents:
- `modis_truecolor`: MODIS true color image
- `modis_falsecolor`: MODIS false color image
- `modis_landmask`: MODIS landmask image
- `modis_cloudfraction`: MODIS cloud fraction image
- `masie_landmask`: MASIE landmask image
- `masie_seaice`: MASIE sea ice image
- `validated_binary_floes`: binary image of validated floes
- `validated_labeled_floes`: labeled image of validated floes
- `validated_floe_properties`: CSV file of validated floe properties

The `dataset` can be iterated over to get each `Case`:
Example:
```jldoctest Watkins2026Dataset
julia> for case in dataset
           println(name(case) * 
                   ": sea ice fraction: " * string(metadata(case).sea_ice_fraction) *
                   ", true color image size: " * string(size(modis_truecolor(case))))
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

```jldoctest Watkins2026Dataset
julia> first(dataset)
Case(GitHubLoader("https://github.com/danielmwatkins/ice_floe_validation_dataset/", "b865acc62f223d6ff14a073a297d682c4c034e5d", "/tmp/Watkins2026"), DataFrameRow
 Row │        case_number  region      start_date  center_lon  center_lat  center_x  center_y  month  sea_ice_fraction  mean_sea_ice_concentration  init_case_number  satellite  visible_sea_ice  visible_la ⋯
     │ Int64  Int64        String      Dates.Date  Float64     Float64     Int64     Int64     Int64  Float64           Float64                     Int64             String     String           String     ⋯
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │    20           11  baffin_bay  2011-07-02    -70.7347     72.3303   -837500  -1737500      7               0.8                       0.388                -1  aqua       yes              no         ⋯
                                                                                                                                                                                            16 columns omitted)

julia> validated_labeled_floes(first(dataset))
Segmented Image with:
labels map: 400×400 Matrix{Int64}
number of labels: 105

julia> modis_truecolor(first(dataset))
400×400 Array{RGBA{N0f8},2} with eltype ColorTypes.RGBA{FixedPointNumbers.N0f8}:
 RGBA{N0f8}(0.094,0.133,0.169,1.0)  RGBA{N0f8}(0.051,0.094,0.118,1.0) ...

```

!!! tip "Cacheing"
    Data are downloaded to the `<cache_dir>/<ref>`, e.g. `/tmp/Watkins2026/b865acc62f223d6ff14a073a297d682c4c034e5d/`. 
    If a file of the correct name already exists in that path, if loaded again the cached data will be returned.

    There are no checks to ensure that the cached data are up-to-date, 
    so if the data change in the source for that `ref`, the loader won't load the new data.
    In this case, you can clear the cache by deleting the cache directory, 
    e.g. `rm -r /tmp/Watkins2026/b865acc62f223d6ff14a073a297d682c4c034e5d/`.

```
"""
function Watkins2026Dataset(;
    url="https://github.com/danielmwatkins/ice_floe_validation_dataset/",
    ref="b865acc62f223d6ff14a073a297d682c4c034e5d",
    cache_dir="/tmp/Watkins2026",
    metadata_path="data/validation_dataset/validation_dataset.csv",
)::Dataset
    loader = GitHubLoader(; url, ref, cache_dir)
    metadata = metadata_path |> loader |> load |> DataFrame
    return Dataset(loader, metadata)
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
    metadata(case).fl_analyst == "" && return nothing
    (; case_number, region, date, satellite) = _filename_parts(case)
    file = "data/validation_dataset/binary_floes/$(case_number)-$(region)-$(date)-$(satellite)-binary_floes.png"
    img = file |> case.loader |> load .|> Gray |> (x -> x .> 0.5) .|> Gray
    return img
end

function validated_labeled_floes(case::Case; ext="tiff")
    metadata(case).fl_analyst == "" && return nothing
    (; case_number, region, date, satellite) = _filename_parts(case)
    file = "data/validation_dataset/labeled_floes/$(case_number)-$(region)-$(date)-$(satellite)-labeled_floes.$(ext)"
    labels = file |> case.loader |> load .|> Int
    img = SegmentedImage(modis_truecolor(case), labels)
    return img
end

function validated_floe_properties(case::Case)::DataFrame
    metadata(case).fl_analyst == "" && return nothing
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

function name(case::Case)::String
    (; case_number, region, date, satellite, pixel_scale, image_scale) = _filename_parts(
        case
    )
    return "$(case_number)-$(region)-$(image_scale)-$(date)-$(satellite)-$(pixel_scale)"
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
