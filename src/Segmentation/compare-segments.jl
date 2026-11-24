"""
    objectwise_comparison(df1, df2, labels1, labels2;
        minimum_area_fraction=0.05,
    )
    objectwise_comparison(labels1, labels2;
        properties=[
            :label, :area, :row_centroid, :col_centroid,
            :max_col, :max_row, :min_col, :min_row
        ],
        minimum_area_fraction=0.05,
    )

Produce a comparison dataframe between connected components in labels1 and labels2. If property tables
have already been computed, all columns will be included, otherwise only columns listed in `properties`
will be included. Column names in df1 will have `1` appended and column names in df2 will have `2`
appended. Note that there may be many objects with no overlap, so the number of rows in the table will
likely not match the number of objects.

The function computes the area overlap and the objectwise precision and recall, where recall is the area 
overlap normalized by `area1` and precision is the overlap normalized by `area2`. All pairs with 
either precision or recall greater than `minimum_area_fraction` are included in the comparison dataframe. 

Three additional similarity metrics are available: `euclidean_distance` (in pixel units), 
`includes_centroid`, which returns Bool indicating if the centroid in one object is included in the other,
and `symmetric_difference`, which calculates the sum of non-overlapping pixels.

Inputs
- `df1`: DataFrame returned from `regionprops_table`. Needs to include at least label, area, bbox, and centroid.
- `df2`: DataFrame returned from `regionprops_table`. Needs to include at least label, area, bbox, and centroid.
- `labels1`: Labeled connected components
- `labels2`: Labeled connected components
- `properties`: List of properties accepted by `regionprops_table` (Only needed if `df1` and `df2` are not supplied)
- `similarity_measures`: List of additional measures of similarity to compute.
- `minimum_area_fraction`: Minimum overlap between objects in `labels1` and `labels2` to include in the table.

Output
- DataFrame with columns for each property from the two input images along with computed similarity metrics.
"""
function objectwise_comparison(
    df1::S, df2::S, labels1::T, labels2::T; 
    similarity_measures::Union{Nothing,Vector{<:AbstractString},Vector{<:Symbol}}=nothing,
    minimum_area_fraction::Real=0.05,
) where {S<:AbstractDataFrame,T<:AbstractMatrix{<:Integer}}

    index1 = component_indices(labels1)
    rs = Dict{Int64,Vector{Int64}}()

    # Collect all instances of overlaps between labels1 and labels2
    for label in keys(index1)
        if label > 0    
            matches = filter(r -> r != 0, unique(labels2[index1[label]]))
            if length(matches) > 0
                push!(rs, label => matches)
            end
        end
    end

    # Add suffixes to distinguish sources
    df1 = rename(n -> n * "1", df1)
    df2 = rename(n -> n * "2", df2)

    # Convert relevant set to DataFrame
    label_lookup = DataFrame([(label1=k, label2=v) for (k, vs) in rs for v in vs])
    df_compare = innerjoin(label_lookup, df1; on=:label1)
    df_compare = innerjoin(df_compare, df2; on=:label2)

    eltype(similarity_measures) <: AbstractString && (similarity_measures = Symbol.(similarity_measures))
    
    transform!(df_compare, AsTable(:) => 
        ByRow(r -> _area_overlap(r; labels1=labels1, labels2=labels2)) =>
            :area_overlap
    )

    transform!(df_compare, [:area1, :area_overlap] => ((A, AUB) -> AUB ./ A) => :recall)
    transform!(df_compare, [:area2, :area_overlap] => ((A, AUB) -> AUB ./ A) => :precision)
    subset!(df_compare, [:recall, :precision] => (r, p) -> max.(r, p) .> minimum_area_fraction)

    isnothing(similarity_measures) && return df_compare

    if :includes_centroid ∈ similarity_measures
        transform!(df_compare, [:row_centroid1, :col_centroid1, :label2] =>
            ByRow((rc, cc, l2) -> _contains_point(rc, cc, l2; labels=labels2)) =>
            :centroid1_in_label2
        )
        transform!(df_compare, [:row_centroid2, :col_centroid2, :label1] =>
            ByRow((rc, cc, l1) -> _contains_point(rc, cc, l1; labels=labels1)) =>
            :centroid2_in_label1
        )
    end

    if :euclidean_distance ∈ similarity_measures
        transform!(df_compare, AsTable(:) =>
            ByRow(r -> _euclidean_distance(r)) =>
            :euclidean_distance
        )
    end


    if :symmetric_difference ∈ similarity_measures
        transform!(df_compare, AsTable(:) =>
            ByRow(r -> _symmetric_difference(r; labels1=labels1, labels2=labels2)) =>
            :symmetric_difference
        )
    end
    
    return df_compare
end

function objectwise_comparison(
    labels1::AbstractMatrix{<:Integer}, labels2::AbstractMatrix{<:Integer}; 
    properties=[:label, :area, :centroid, :bbox],
    similarity_measures::Union{Nothing,Vector{<:AbstractString},Vector{<:Symbol}}=nothing,
    minimum_area_fraction::Real=0.05
    )
    df1 = regionprops_table(labels1; properties=properties)
    df2 = regionprops_table(labels2; properties=properties)
    return objectwise_comparison(df1, df2, labels1, labels2; 
        similarity_measures=similarity_measures,
        minimum_area_fraction=minimum_area_fraction,
        )
end

### Similarity Measures
# Each expects the input to be either a named tuple or a DataFrameRow
# with columns `a1`, `a2`, `b1`, `b2` where `a`, `b` are properties returned
# by the `regionprops_table`. Some similarity measures may need additional inputs.

"""
    _area_overlap(df_row; labels1, labels2)

Helper function to compute the area overlap between labeled regions using a row in the comparison dataframe.
"""
function _area_overlap(df_row; 
    labels1::AbstractMatrix{<:Integer}, 
    labels2::AbstractMatrix{<:Integer})
    rmin = min(df_row.min_row1, df_row.min_row2)
    rmax = max(df_row.max_row1, df_row.max_row2)
    cmin = min(df_row.min_col1, df_row.min_col2)
    cmax = max(df_row.max_col1, df_row.max_col2)
    mask1 = labels1[rmin:rmax, cmin:cmax] .== df_row.label1
    mask2 = labels2[rmin:rmax, cmin:cmax] .== df_row.label2
    intersect_area = sum(mask1 .&& mask2)
    return intersect_area
end

"""
    _contains_point(row, col, label; labels::Matrix{Integer})

Helper function to determine whether a point (row, col) is contained within an object labeled `label`
in the indexmap `labels`. 

"""
function _contains_point(
    row::Real, col::Real, label::Integer; 
    labels::AbstractMatrix{<:Integer})::Bool
    row = round(Int64, row)
    col = round(Int64, col)
    return labels[row, col] == label
end

"""
    _euclidean_distance(df_row)

Helper function to compute euclidean distance in pixels for a comparison dataframe.
"""
function _euclidean_distance(df_row)
    r1 = df_row.row_centroid1
    r2 = df_row.row_centroid2
    c1 = df_row.col_centroid1
    c2 = df_row.col_centroid2
    return sqrt((r1 - r2)^2 + (c1 - c2)^2)
end


"""
    _symmetric_difference(df_row; labels1, labels2)

Helper function to compute the symmetric difference between labeled regions
using a row in the comparison dataframe.
"""
function _symmetric_difference(df_row; 
    labels1::AbstractMatrix{<:Integer}, 
    labels2::AbstractMatrix{<:Integer})
    rmin = min(df_row.min_row1, df_row.min_row2)
    rmax = max(df_row.max_row1, df_row.max_row2)
    cmin = min(df_row.min_col1, df_row.min_col2)
    cmax = max(df_row.max_col1, df_row.max_col2)
    mask1 = labels1[rmin:rmax, cmin:cmax] .== df_row.label1
    mask2 = labels2[rmin:rmax, cmin:cmax] .== df_row.label2
    return sum(mask1 .!= mask2)
end


### Relevant set can be computed from the comparison dataframe
# TODO: Add option to return the computed area overlap to prevent repeated computations, e.g. as DataFrame
"""
    get_relevant_set(df1, df2, labels1, labels2; minimum_area_fraction=0.5)
    get_relevant_set(labels1, labels2; minimum_area_fraction=0.5)

Find the relevant set for comparing two segmentation results. The relevant set for a
segmentation comparison set s in S in reference to object g in G is defined by

1. centroid g in s
2. centroid s in g
3. area overlap greater than 50% of g
4. area overlap greater than 50% of s

Returns a DataFrame with columns `label1`, `label2`.

Inputs:
- `df1`, `df2` = results of regionprops table (optional)
- `labels1`, `labels2` = image indexmaps
- `minimum_area_fraction` = 0.5 Minimum overlap fraction to count as part of the relevant set

References:
Clinton et al. (2010) "Accuracy Assessment Measures for Object-based Image Segmentation Goodness," Photogrammetric Engineering and Remote Sensing, 76, 3, pp 289-299. [Link](https://www.ingentaconnect.com/content/asprs/pers/2010/00000076/00000003/art00004)

"""
function get_relevant_set(df1, df2, labels1, labels2; minimum_area_fraction=0.5)::DataFrame

    # First get all the comparison points
    df_comparison = objectwise_comparison(
        df1, df2, labels1, labels2; 
            similarity_measures=[:includes_centroid], 
            minimum_area_fraction=0
    )

    transform!(df_comparison, [:precision, :recall] =>
        ByRow((p, r) -> max.(p, r) > minimum_area_fraction) =>
        :area_fraction_test
    )
    
    # Then check the centroids
    transform!(df_comparison, [:centroid1_in_label2, :centroid2_in_label1] =>
        ByRow((c1, c2) -> max.(c1, c2)) => :centroid_test
    )

    # If either test passes, it is part of the relevant set
    subset!(df_comparison, [:area_fraction_test, :centroid_test] =>
        ((A, C) -> A .|| C))

    # All remaining objects are part of the relevant sets
    return df_comparison[:, [:label1, :label2]] 
end

function get_relevant_set(labels1, labels2; minimum_area_fraction=0.5)::DataFrame
    df1 = regionprops_table(labels1; properties=[:label, :area, :bbox, :centroid])
    df2 = regionprops_table(labels2; properties=[:label, :area, :bbox, :centroid])
    return get_relevant_set(df1, df2, labels1, labels2; minimum_area_fraction=minimum_area_fraction)
end
