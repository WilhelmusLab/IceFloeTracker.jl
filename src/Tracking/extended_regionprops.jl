# Functions for adding additional columns to regionprops needed for floe tracking

import DataFrames: DataFrame, nrow
import Images: label_components
import ..Morphology: bwareamaxfilt

FloeLabelsImage = Union{BitMatrix, Matrix{<:Bool}, Matrix{<:Integer}}

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
                props_row.label
            )

    elseif issubset(bbox_column_names, colnames)
        floesimg_bitmatrix = floesimg .> 0
        return cropfloe(
            floesimg_bitmatrix,
            props_row.min_row,
            props_row.min_col,
            props_row.max_row,
            props_row.max_col
        )
    
    elseif issubset(label_column_names, colnames)
        return cropfloe(floesimg, props_row.label)
    
    end
end

"""
    cropfloe(floesimg, min_row, min_col, max_row, max_col)

Crops the floe delimited by `min_row`, `min_col`, `max_row`, `max_col`, from the floe image `floesimg`.
"""
function cropfloe(floesimg::BitMatrix, min_row::I, min_col::I, max_row::I, max_col::I) where {I<:Integer}
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
function cropfloe(floesimg::Matrix{I}, min_row::J, min_col::J, max_row::J, max_col::J, label::I)  where {I<:Integer, J<:Integer}
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
    addfloemasks!(props::DataFrame, floeimg::FloeLabelsImage)

Add a column to `props` called `floearray` containing the cropped floe masks from `floeimg`.
"""
function addfloemasks!(props::DataFrame, floeimg::FloeLabelsImage)
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
