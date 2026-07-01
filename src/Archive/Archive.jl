
module Archive

using NCDatasets, Images, Dates, TimeZones, DataFrames
import ..Geospatial: latlon
import ..Segmentation: regionprops_table, converttounits!
import ..ImageUtils: binarize_mask

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

Convert missing values in Float64 columns of the DataFrame `df` to `NaN` to allow saving as HDF5/NetCDF.
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
        VersionNumber(file.attrib["file_version"])
    end
    if version == VersionNumber("3.0.0")
        return _load_v3(input_path)
    else
        error("Unsupported file version: $version")
    end
end

include("./V3.jl")

end
