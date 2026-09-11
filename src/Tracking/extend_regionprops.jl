# Functions for adding additional columns to regionprops needed for floe tracking

import DataFrames: DataFrame, nrow, DataFrameRow, transform!, ByRow, AbstractDataFrame
import Images: label_components, SegmentedImage, labels_map
import Dates: DateTime
import Random: randstring
import ..Segmentation: component_floes

FloeLabelsImage = Union{BitMatrix,Matrix{<:Bool},Matrix{<:Integer},<:SegmentedImage}

# TODO: Change "passtimes" to "image_time". In principle we could be using images from airplanes/helicopters, not just satellites
"""
    add_passtimes!(props::DataFrame, passtimes::DateTime)
    add_passtimes!.(props::Vector{DataFrame}, passtimes::Vector{DateTime})

Add a column `passtime` to each DataFrame in `props` containing the time of the image in which the floes were captured.

# Arguments
- `props`: array of DataFrames containing floe properties.
- `passtimes`: array of `DateTime` objects containing the time of the image in which the floes were captured.

"""
function add_passtimes!(props_df::DataFrame, passtime::DateTime)
    props_df.passtime .= passtime
    return nothing
end

"""
    add_ψs!(props_df::DataFrame})
    add_ψs!.(props_dfs::Vector{DataFrame})

Add the ψ-s curves to each row of `props_df`.

Note: each member of `props` must have a `boundary` column with the traced, resampled
floe boundary. To add boundaries see [`add_boundary!`](@ref), which must run first.
"""
function add_ψs!(props_df::DataFrame)
    hasproperty(props_df, :boundary) || throw(
        ArgumentError("`add_ψs!` requires a `:boundary` column; call `add_boundary!` first"),
    )
    props_df.psi = map(boundary -> buildψs(boundary)[1], props_df.boundary)
    return nothing
end

_uuid() = randstring(12)

"""
    add_uuids!(df::DataFrame)
    add_uuids!.(dfs::Vector{DataFrame})

Assign a unique ID to each floe in a (vector of) table(s) of floe properties.
"""
function add_uuids!(df::DataFrame)
    df.uuid = [_uuid() for _ in 1:nrow(df)]
    return df
end

"""
    add_floemasks!(props::DataFrame, indexmap::FloeLabelsImage)
    add_floemasks!.(props::Vector{DataFrame}, indexmap::Vector{FloeLabelsImage})

Add a column to `props` called `mask` containing the cropped floe masks from `indexmap`.
"""
function add_floemasks!(
    props::DataFrame, indexmap::Matrix{Int64}; label_column::Symbol=:label
)
    floes = component_floes(indexmap)
    img_labels = props[:, label_column]
    props[:, :mask] = map(s -> floes[s], img_labels)
    return nothing
end

function add_floemasks!(
    props::DataFrame, segmented_image::SegmentedImage; label_column::Symbol=:label
)
    add_floemasks!(props, labels_map(segmented_image); label_column=label_column)
    return nothing
end

"""
    add_boundary!(props_df::DataFrame; reduc_factor::Int64=2)
    add_boundary!.(props_dfs::Vector{DataFrame}; reduc_factor::Int64=2)

Add resampled boundary curves as `:boundary` column to DataFrame.
Each boundary is a Matrix{Float64}(n, 2) with [x y] coordinates.

Note: each member of `props_df` must have a `:mask` column with a binary image representing the floe.
To add floe masks see [`add_floemasks!`](@ref).

# Arguments
- `props_df`: DataFrame with `:mask` column
- `reduc_factor`: Reduction factor for boundary resampling (default: 2 = 50% reduction)
"""
function add_boundary!(props_df::DataFrame; reduc_factor::Int64=2)
    props_df.boundary = map(mask -> _traced_boundary(mask; reduc_factor), props_df.mask)
    return nothing
end

"""
    _traced_boundary(mask; reduc_factor::Int64=2)

Trace a floe boundary and resample it by arc length.

Shared by [`add_boundary!`](@ref) and `buildψs(::AbstractArray)` so the two cannot drift
apart.
"""
function _traced_boundary(mask; reduc_factor::Int64=2)
    bd_traced = bwtraceboundary(mask)
    # Handle case of multiple boundaries (shouldn't happen for individual floes but be defensive)
    bd_traced_single = isa(bd_traced, Vector{Vector{CartesianIndex}}) ?
        bd_traced[1] : bd_traced
    return resample_boundary(bd_traced_single, reduc_factor)
end
