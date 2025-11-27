# Functions for adding additional columns to regionprops needed for floe tracking

import DataFrames: DataFrame, nrow, DataFrameRow, transform!, ByRow, AbstractDataFrame
import Images: label_components, SegmentedImage, labels_map
import Dates: DateTime
import ..Morphology: bwareamaxfilt

FloeLabelsImage = Union{BitMatrix,Matrix{<:Bool},Matrix{<:Integer},<:SegmentedImage}

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

Note: each member of `props` must have a `mask` column with a binary image representing the floe. 
To add floe masks see [`addfloemasks!`](@ref).
"""
function add_ψs!(props_df::DataFrame)
    props_df.psi = map(buildψs, props_df.mask)
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

# TODO: Update the cropfloes function to use the "label" parameter in the regionprops table.
# This way, we can create a bitmatrix with labeled image == label, and crop that.
# TODO: Add method to allow SegmentedImage as input
# TODO: bbox and label names as keyword arguments
"""
    cropfloe(floesimg, props, i)

Crops the floe delimited by the bounding box data in `props` at index `i` from the floe image `floesimg`.

If the dataframe has bounding box data `min_row`, `min_col`, `max_row`, `max_col`, but no `label`, then returns the largest contiguous component.

If the dataframe has bounding box data `min_row`, `min_col`, `max_row`, `max_col`, and a `label`, then returns the component with the label. In this case, `floesimg` must be an Array{Int}.

If the dataframe has only a `label` and no bounding box data, then returns the component with the label, padded by one cell of zeroes on all sides. In this case, `floesimg` must be an Array{Int}.


"""
function cropfloe(floesimg::FloeLabelsImage, props::DataFrame, i::Integer)
    props_row = props[i, :]
    colnames = Set(names(props_row))
    bbox_column_names = Set(["min_row", "min_col", "max_row", "max_col"])
    label_column_names = Set(["label"])
    bbox_label_column_names = union(bbox_column_names, label_column_names)

    if issubset(bbox_label_column_names, colnames)
        return cropfloe(
            floesimg,
            props_row.min_row,
            props_row.min_col,
            props_row.max_row,
            props_row.max_col,
            props_row.label,
        )

    elseif issubset(bbox_column_names, colnames)
        floesimg_bitmatrix = floesimg .> 0
        return cropfloe(
            floesimg_bitmatrix,
            props_row.min_row,
            props_row.min_col,
            props_row.max_row,
            props_row.max_col,
        )

    elseif issubset(label_column_names, colnames)
        return cropfloe(floesimg, props_row.label)
    end
end

"""
    cropfloe(floesimg, min_row, min_col, max_row, max_col)

Crops the floe delimited by `min_row`, `min_col`, `max_row`, `max_col`, from the floe image `floesimg`.
"""
function cropfloe(
    floesimg::BitMatrix, min_row::I, min_col::I, max_row::I, max_col::I
) where {I<:Integer}
    #= 
    Crop the floe using bounding box data in props.
    Note: Using a view of the cropped floe was considered but if there were multiple components in the cropped floe, the source array with the floes would be modified. =#
    prefloe = floesimg[min_row:max_row, min_col:max_col]

    #= Check if more than one component is present in the cropped image.
    If so, keep only the largest component by removing all on pixels not in the largest component =#
    components = label_components(prefloe, trues(3, 3))

    if length(unique(components)) > 2
        mask = bwareamaxfilt(components .> 0)
        prefloe[.!mask] .= 0
    end
    return prefloe
end

"""
    cropfloe(floesimg, min_row, min_col, max_row, max_col, label)

Crops the floe from `floesimg` with the label `label`, returning the region bounded by `min_row`, `min_col`, `max_row`, `max_col`, and converting to a BitMatrix.
"""
function cropfloe(
    floesimg::Matrix{I}, min_row::J, min_col::J, max_row::J, max_col::J, label::I
) where {I<:Integer,J<:Integer}
    #= 
    Crop the floe using bounding box data in props.
    Note: Using a view of the cropped floe was considered but if there were multiple components in the cropped floe, the source array with the floes would be modified. =#
    prefloe = floesimg[min_row:max_row, min_col:max_col]
    @debug "prefloe: $prefloe"

    #= Remove any pixels not corresponding to that numbered floe 
    (each segment has a different integer) =#
    floe_area = prefloe .== label
    @debug "mask: $floe_area"

    return floe_area
end

"""
    add_floemasks!(props::DataFrame, floeimg::FloeLabelsImage)
    add_floemasks!.(props::Vector{DataFrame}, floeimgs::Vector{FloeLabelsImage})

Add a column to `props` called `mask` containing the cropped floe masks from `floeimg`.
"""
function add_floemasks!(props::DataFrame, floeimg::FloeLabelsImage)
    props.mask = getfloemasks(props, floeimg)
    return nothing
end

"""
    getfloemasks(props::DataFrame, floeimg::BitMatrix)

Return a vector of cropped floe masks from `floeimg` using the bounding box data in `props`.
"""
function getfloemasks(props::DataFrame, floeimg::FloeLabelsImage)
    return map(i -> cropfloe(floeimg, props, i), 1:nrow(props))
end

function add_floemasks!(
    props::DataFrame, segmented_image::SegmentedImage; label_column::Symbol=:label
)
    floeimg = labels_map(segmented_image)
    props[!, :mask] = missings(BitMatrix, nrow(props))

    for i in 1:nrow(props)
        label = props[i, label_column]
        min_row = props[i, :min_row]
        min_col = props[i, :min_col]
        max_row = props[i, :max_row]
        max_col = props[i, :max_col]
        props.mask[i] = cropfloe(floeimg, min_row, min_col, max_row, max_col, label)
    end
    return props
end
