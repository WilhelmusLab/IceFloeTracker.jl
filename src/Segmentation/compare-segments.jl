"""
    objectwise_comparison(df1, df2, labels1, labels2;
        tol_area_fraction=0.05,
    )
    objectwise_comparison(labels1, labels2;
        properties=[
            :label, :area, :row_centroid, :col_centroid,
            :max_col, :max_row, :min_col, :min_row
        ],
        tol_area_fraction=0.05,
    )

Produce a comparison dataframe between connected components in labels1 and labels2. If property tables
have already been computed, all columns will be included, otherwise only columns listed in `properties`
will be included. All pairs with area overlap greater than `tol_area_fraction` are included in the
comparison dataframe. Column names in df1 will have `1` appended and column names in df2 will have `2`
appended. Note that there may be many objects with no overlap, so the number of rows in the table will
likely not match the number of objects.

In addition to the alignment of the standard object properties, the function computes 4 additional measures
of similarity: the Euclidean distance between centroids; object-wise precision and recall; and the area of
overlap between the two labels.

Inputs
- `df1`: DataFrame returned from `regionprops_table`
- `df2`: DataFrame returned from `regionprops_table`
- `labels1`: Labeled connected components
- `labels2`: Labeled connected components
- `properties`: List of properties accepted by `regionprops_table`
- `tol_area_fraction`: Minimum overlap between objects in `labels1` and `labels2` to include in the table.

Output
- DataFrame with columns for each property from the two input images along with computed similarity metrics.
"""
function objectwise_comparison(
    df1::S, df2::S, labels1::T, labels2::T; tol_area_fraction::Real=0.05
) where {S<:AbstractDataFrame,T<:AbstractMatrix{<:Integer}}

    rs = get_relevant_set(df1, df2, labels1, labels2; minimum_area_fraction=tol_area_fraction)

    # Add suffixes to distinguish sources
    df1 = rename(n -> n * "1", df1)
    df2 = rename(n -> n * "2", df2)

    # Convert relevant set to DataFrame
    label_lookup = DataFrame([(label1=k, label2=v) for (k, vs) in rs for v in vs])
    df_compare = innerjoin(label_lookup, df1; on=:label1)
    df_compare = innerjoin(df_compare, df2; on=:label2)

    transform!(df_compare,
        [:row_centroid1, :row_centroid2, :col_centroid1, :col_centroid2] =>
            ByRow((r1, r2, c1, c2) -> sqrt((r1 - r2)^2 + (c1 - c2)^2)) =>
                :euclidean_distance
    )

    # AsTable passes the dataframe as a named tuple, then ByRow iterates through the tuple.
    transform!(df_compare, AsTable(:) => 
        ByRow(r -> _area_overlap(r; labels1=labels1, labels2=labels2)) =>
            :area_overlap
    )

    transform!(df_compare, [:area1, :area_overlap] => ByRow((A, AUB) -> AUB/A) => :recall)
    transform!(df_compare, [:area2, :area_overlap] => ByRow((A, AUB) -> AUB/A) => :precision)

    return df_compare
end

function objectwise_comparison(labels1, labels2; 
    properties=[
            :label, :area, :centroid, :bbox
    ],
    tol_area_fraction=0.05
    )
    df1 = regionprops_table(labels1; properties=properties)
    df2 = regionprops_table(labels2; properties=properties)
    return objectwise_comparison(df1, df2, labels1, labels2; tol_area_fraction=tol_area_fraction)
end

"""
    _area_overlap(df_row; labels1, labels2)

Helper function to compute the area overlap between labeled regions using a row in the comparison dataframe.
"""
function _area_overlap(df_row; labels1::Matrix{Int64}, labels2::Matrix{Int64})
    rmin = min(df_row.min_row1, df_row.min_row2)
    rmax = max(df_row.max_row1, df_row.max_row2)
    cmin = min(df_row.min_col1, df_row.min_col2)
    cmax = max(df_row.max_col1, df_row.max_col2)
    mask1 = labels1[rmin:rmax, cmin:cmax] .== df_row.label1
    mask2 = labels2[rmin:rmax, cmin:cmax] .== df_row.label2
    intersect_area = sum(mask1 .&& mask2)
    return intersect_area
end

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

Returns a dictionary mapping labels in `labels1` to a set of labels in `labels2`.

Inputs:
- `df1`, `df2` = results of regionprops table (optional)
- `labels1`, `labels2` = image indexmaps
- `minimum_area_fraction` = 0.5 Minimum overlap fraction to count as part of the relevant set

References:
Clinton et al. (2010) "Accuracy Assessment Measures for Object-based Image Segmentation Goodness," Photogrammetric Engineering and Remote Sensing, 76, 3, pp 289-299. [Link](https://www.ingentaconnect.com/content/asprs/pers/2010/00000076/00000003/art00004)

"""
function get_relevant_set(df1, df2, labels1, labels2; minimum_area_fraction=0.5)::Dict{Int64,Vector{Int64}}
    relevant_set = Dict{Int64,Vector{Int64}}()
    for floe in eachrow(df1)
        # select labels that are inside the bounding box for the floe
        matched_labels = unique(
            labels2[floe.min_row:floe.max_row, floe.min_col:floe.max_col]
        )

        # if any, then check centroid positions
        maximum(matched_labels) == 0 && continue
        # get the rows in the segments_df from the matched labels
        candidate_subset = subset(df2, :label => ByRow(l -> l in matched_labels))

        relevant_set_labels = Int64[]

        # check if centroid g in s
        rc = round(Int64, floe.row_centroid)
        cc = round(Int64, floe.col_centroid)
        push!(relevant_set_labels, labels2[rc, cc])

        # check if centroid s in g
        for s_floe in eachrow(candidate_subset)
            rc = round(Int64, s_floe.row_centroid)
            cc = round(Int64, s_floe.col_centroid)
            (labels1[rc, cc] == floe.label) && begin
                push!(relevant_set_labels, s_floe.label)
            end

            # joint bbox
            rmin = min(floe.min_row, s_floe.min_row)
            rmax = max(floe.max_row, s_floe.max_row)
            cmin = min(floe.min_col, s_floe.min_col)
            cmax = max(floe.max_col, s_floe.max_col)

            # check if area overlap between g and s is larger than 50% of g
            gtmask = labels1[rmin:rmax, cmin:cmax] .== floe.label
            slmask = labels2[rmin:rmax, cmin:cmax] .== s_floe.label
            intersect_area = sum(gtmask .&& slmask)
            if max(intersect_area / s_floe.area, intersect_area / floe.area) > minimum_area_fraction
                push!(relevant_set_labels, s_floe.label)
            end
        end
        relevant_set_labels = filter(r -> r != 0, unique(relevant_set_labels))
        if length(relevant_set_labels) > 0
            push!(relevant_set, floe.label => relevant_set_labels)
        end
    end
    return relevant_set
end

function get_relevant_set(labels1, labels2; minimum_area_fraction=0.5)::Dict{Int64,Vector{Int64}}
    df1 = regionprops_table(labels1; properties=[:label, :area, :bbox, :centroid])
    df2 = regionprops_table(labels2; properties=[:label, :area, :bbox, :centroid])
    return get_relevant_set(df1, df2, labels1, labels2; minimum_area_fraction=minimum_area_fraction)
end
