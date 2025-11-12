# Functions for adding additional columns to regionprops needed for floe tracking

import DataFrames: DataFrame, nrow, DataFrameRow, transform!, ByRow, AbstractDataFrame
import Images: label_components
import ..Morphology: bwareamaxfilt

FloeLabelsImage = Union{BitMatrix, Matrix{<:Bool}, Matrix{<:Integer}}
abstract type AbstractThresholdFunction <: Function end


"""
    add_passtimes!(props, passtimes)

Add a column `passtime` to each DataFrame in `props` containing the time of the image in which the floes were captured.

# Arguments
- `props`: array of DataFrames containing floe properties.
- `passtimes`: array of `DateTime` objects containing the time of the image in which the floes were captured.

"""
function add_passtimes!(props, passtimes)
    for (i, passtime) in enumerate(passtimes)
        props[i].passtime .= passtime
    end
    return nothing
end


"""
    addψs!(props::Vector{DataFrame})

Add the ψ-s curves to each member of `props`.

Note: each member of `props` must have a `mask` column with a binary image representing the floe. 
To add floe masks see [`addfloemasks!`](@ref).
"""
function addψs!(props::Vector{DataFrame})
    for prop in props
        prop.psi = map(buildψs, prop.mask)
    end
    return nothing
end

"""
    addψs!(props_df::DataFrame})

Add the ψ-s curves to each row of `props_df`.

Note: each member of `props` must have a `mask` column with a binary image representing the floe. 
To add floe masks see [`addfloemasks!`](@ref).
"""
function addψs!(props_df::DataFrame)
    props_df.psi = map(buildψs, props_df.mask)
    return nothing
end

function addfloemasks!(props::Vector{DataFrame}, imgs::Vector{<:FloeLabelsImage})
    for (img, prop) in zip(imgs, props)
        addfloemasks!(prop, img)
    end
    return nothing
end

_uuid() = randstring(12)

"""
    adduuid!(df::DataFrame)
    adduuid!(dfs::Vector{DataFrame})

Assign a unique ID to each floe in a (vector of) table(s) of floe properties.
"""
function adduuid!(df::DataFrame)
    df.uuid = [_uuid() for _ in 1:nrow(df)]
    return df
end

function adduuid!(dfs::Vector{DataFrame})
    for (i, _) in enumerate(dfs)
        adduuid!(dfs[i])
    end
    return dfs
end

"""
    _add_integer_id!(df, col, new)

For distinct values in the column `col` of `df`, add a new column `new` to be consecutive integers starting from 1.
"""
function _add_integer_id!(df::AbstractDataFrame, col::Symbol, new::Symbol)
    ids = unique(df[!, col])
    _map = Dict(ids .=> 1:length(ids))
    transform!(df, col => ByRow(x -> _map[x]) => new)
    return nothing
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
""" #TODO: have the names of the columns created in the function be named parameters
function time_distance_test!(
        floe::DataFrameRow,
        candidates::DataFrame;
        threshold_function=LopezAcostaTimeDistanceFunction(),
        threshold_column=:time_distance_test
        )
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
    sd(mask, orientation) = round(shape_difference(floe.mask, floe.orientation, mask, orientation), digits=3)

    transform!(candidates,  [:mask, :orientation] => 
        ByRow(sd) => :shape_difference
        )

    candidates[!, :scaled_shape_difference] = candidates[!, :shape_difference] ./ candidates[!, scale_by]
    candidates[!, :scaled_shape_difference] .= round.(candidates[!, :scaled_shape_difference], digits=3)

    transform!(candidates, [area_column, :scaled_shape_difference] =>
        ByRow(threshold_function) => threshold_column
    )
end

#TODO: Add option to include the confidence intervals with the normalized cross correlation tests.

"""
    psi_s_correlation_test!(floe, candidates;
     threshold_function, threshold_column=:psi_s_correlation_test, area_column=:area)

   Compute the psi-s correlation between a floe and a dataframe of candidate floes. Adds the 
   psi-s correlation,  psi-s correlation score (1 - correlation), and the result of the threshold function
   to the columns of `candidates`.
""" 
function psi_s_correlation_test!(
    floe::DataFrameRow,
    candidates::DataFrame;
    threshold_function::Function,
    threshold_column=:psi_s_correlation_test,
    area_column=:area
)
    if :psi ∉ names(candidates)
        p1 = buildψs(floe.mask)
        addψs!(candidates)
    else
        p1 = floe.psi
    end
    rfloe(p2) = round(normalized_cross_correlation(p1, p2), digits=3)
    transform!(candidates,  [:psi] => ByRow(rfloe) => :psi_s_correlation)
    candidates[!, :psi_s_correlation_score] = 1 .- candidates[!, :psi_s_correlation]

    # Future work: add computation of the confidence intervals for psi-s corr here.

    transform!(candidates, [area_column, :psi_s_correlation_score] =>
        ByRow(threshold_function) => threshold_column
    )
    subset!(candidates, Not(:psi))
end