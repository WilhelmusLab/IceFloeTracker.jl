# Functions for adding additional columns to regionprops needed for floe tracking

import DataFrames: DataFrame, nrow, DataFrameRow, transform!, ByRow
import Images: label_components
import ..Morphology: bwareamaxfilt

FloeLabelsImage = Union{BitMatrix, Matrix{<:Bool}, Matrix{<:Integer}}

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




####### Extend + Threshold Test Functions #######
# These functions add columns to the candidate matches table and
# perform tests against threshold functions.

abstract type AbstractThresholdFunction end

# This new function also works for the LopezAcosta2019 functions
# as long as you allow for some small amount in between the maximum and minimum area.
"""
The piecewise linear threshold function is defined using two (area, value) pairs. For
areas below the minimum area, it is constant at minimum value; likewise for above the
maximum area. The threshold function is linear in between these two points. A return 
value `true` indicates that the value is below the threshold. 
"""
@kwdef struct PiecewiseLinearThresholdFunction <: AbstractThresholdFunction
    minimum_area = 100
    maximum_area = 700
    minimum_value = 0.4
    maximum_value = 0.2
end

function (f::PiecewiseLinearThresholdFunction)(area, value)
    area < f.minimum_area && return value < f.minimum_value
    area > f.maximum_area && return value < f.maximum_value
    slope = (f.maximum_value - f.minimum_value) / (f.maximum_area - f.minimum_area)
    return value < slope*(area - f.maximum_area) + f.maximum_value
end

"""
    relative_error_test!(floe::DataFrameRow,
     candidates::DataFrame,
     variable::Union{Symbol,String},
     new_column::Union{Symbol,String},
     area_variable::Union{Symbol,String},
     theshold_function::Function
    )

Compute and test (absolute) relative error for `variable`. The relative error
between scalar variables X and Y is defined as 
```
err = abs(X - Y)/mean(X, Y)
```
This function takes a scalar `<variable>` (which must be a named column in 
the `candidates` DataFrame) and computes the relative error. It adds a new
column with the pattern `relative_error_<variable>` and uses a threshold function
to determine if the relative error is small enough given the floe area. Mutates
the candidate dataframe in place.

"""
function relative_error_test!(
     floe::DataFrameRow,
     candidates::DataFrame;
     variable::Union{Symbol,String},
     threshold_column::Union{Symbol,String},
     area_variable::Union{Symbol,String},
     threshold_function::Union{Function, AbstractThresholdFunction}
    )
    new_variable = Symbol(:relative_error_, variable)
    X = floe[variable]
    Y = candidates[!, variable]
    candidates[!, new_variable] = abs.(X .- Y) ./ (0.5 .* (X .+ Y))
    transform!(candidates, [area_variable, new_variable] => ByRow(threshold_function) => threshold_column)
end

"""
    euclidean_distance(x, Y; r)

Compute Euclidean distance between a floe and a dataframe of candidate floes assuming pixels of size `r`.
Each column needs to have centroids calculated already.
"""
function euclidean_distance(floe, candidates; r = 250)
    return sqrt.((floe.row_centroid .- candidates.row_centroid).^2 .+ 
    (floe.col_centroid .- candidates.col_centroid).^2) * r
end

"""
    time_distance_test(floe, candidates; threshold_function, threshold_column)
"""
function time_distance_test!(
        floe::DataFrameRow,
        candidates::DataFrame;
        threshold_function=LopezAcostaTimeDistanceFunction(),
        threshold_column=:time_distance_test)
    candidates[!, :Δt] = candidates[!, :passtime] .- floe.passtime
    candidates[!, :Δx] = euclidean_distance(floe, candidates)
    transform!(candidates, [:Δx, :Δt] => 
        ByRow(threshold_function) => threshold_column)
end

"""
    shape_difference_test(floe, candidates;
     threshold_function, threshold_column=:shape_difference_test, scale_by=:area, area_column=:area)

    Test the scaled shape difference between input `floe` and each floe in the dataframe `candidates`.
    Assumes that the shape difference test operates on the shape difference scaled by a variable `scale_by`
    and the shape difference test depends on the area.

""" 
function shape_difference_test!(
    floe::DataFrameRow,
    candidates::DataFrame;
    threshold_function::Function,
    threshold_column=:shape_difference_test,
    scale_by=:area,
    area_column=:area
)
    sd(mask, orientation) = shape_difference(floe.mask, floe.orientation, mask, orientation)

    transform!(candidates,  [:mask, :orientation] => 
        ByRow(sd) => :shape_difference
        )

    candidates[!, :scaled_shape_difference] = candidates[!, :shape_difference] ./ candidates[!, scale_by]

    transform!(candidates, [area_column, :scaled_shape_difference] =>
        ByRow(threshold_function) => threshold_column
    )
end


"""
    psi_s_correlation_test!(floe, candidates;
     threshold_function, threshold_column=:psi_s_correlation_test, area_column=:area)

    Compute the psi-s correlation between a floe and a dataframe of candidate floes. 

""" 
function psi_s_correlation_test!(
    floe::DataFrameRow,
    candidates::DataFrame;
    threshold_function::Function,
    threshold_column=:psi_s_correlation_test,
    area_column=:area
)
    addψs!(candidates)

    candidates[!, :psi_s_correlation] 
    candidates[!, :psi_s_correlation_score] = 1 .- candidates[!, :psi_s_correlation]

    transform!(candidates, [area_column, :psi_s_correlation_score] =>
        ByRow(threshold_function) => threshold_column
    )
end