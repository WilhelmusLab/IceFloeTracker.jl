"""
    Archive

Module for saving and loading IceFloeTracker.jl segmentation results as netCDF-4
files. The public API is:

- [`save`](@ref) — write a format object (e.g. [`V1`](@ref)) to disk.
- [`load`](@ref) — read a file and return the appropriate format object; the
  format version is detected automatically from the `ift_archive_version` global
  attribute.

Each file-format version is implemented in its own submodule (`ArchiveV1`,
`ArchiveV2`, …) which extends `save` and `load` for its own type. The public
format types (e.g. `V1`) are re-exported from this module so callers only need
`IceFloeTracker.Archive`.
"""
module Archive

using NCDatasets, Images, DataFrames

"""
    save(output_path::AbstractString, obj)

Write an archive object to `output_path` as a netCDF-4 file. Extended by each
archive-format submodule (e.g. `ArchiveV1`) for its own object type.
"""
function save end

"""
    load(input_path::AbstractString, type)

Load an archive object from `input_path` as a netCDF-4 file. Extended by each
archive-format submodule (e.g. `ArchiveV1`) for its own object type.
"""
function load end

function choose_dtype(mx::T) where {T<:Integer}
    types = [UInt8, Int8, UInt16, Int16, UInt32, Int32, UInt64, Int64]
    for t_ in types
        if typemin(t_) <= mx <= typemax(t_)
            return t_
        end
    end
    return error("$mx cannot be represented by any of $types")
end

"""
    convert_missing_to_nan!(df)
    convert_missing_to_nan(df)

Convert missing values in Float64 columns of the DataFrame `df` to `NaN` to allow saving as NetCDF.
"""
function convert_missing_to_nan!(df::DataFrame)
    for (col_name, col_data) in pairs(eachcol(df))
        if eltype(col_data) <: Union{Missing,Float64}
            col_data .= coalesce.(col_data, NaN)
            disallowmissing!(df, col_name)
        end
    end
end

function convert_missing_to_nan(df::DataFrame)
    df_copy = copy(df)
    convert_missing_to_nan!(df_copy)
    return df_copy
end

"""
    load(input_path::AbstractString)

Load an IceFloeTracker.jl archive file. 
"""
function load(input_path::AbstractString)
    version = NCDataset(input_path, "r") do file
        VersionNumber(file.attrib["ift_archive_version"])
    end
    if version == VersionNumber("1.0.0")
        return load(input_path, V1)
    else
        error("Unsupported file version: $version")
    end
end

include("./ArchiveV1.jl")

using .ArchiveV1: V1

end
