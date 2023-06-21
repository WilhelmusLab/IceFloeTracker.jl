
"""
    regionprops_table(label_img, intensity_img; properties, connectivity, extra_properties)

A wrapper of the `regionprops_table` function from the skimage python library.
    
See its full documentation at https://scikit-image.org/docs/stable/api/skimage.measure.html#regionprops-table.
    
# Arguments
- `label_img`: Image with the labeled objects of interest
- `intensity_img`: (Optional) Used for generating `extra_properties`, integer/float array from which (presumably) `label_img` was generated 
- `properties`: List (`Vector` or `Tuple`) of properties to be generated for each connected component in `label_img`
- `extra_properties`: (Optional) not yet implemented. It will be set to `nothing`

# Notes
- Zero indexing has been corrected for the `bbox` and `centroid` properties
- `bbox` data (`max_col` and `max_row`) are inclusive
- `centroid` data are rounded to the nearest integer

See also [`regionprops`](@ref)

# Examples

```jldoctest; setup = :(using IceFloeTracker, Random)
julia> using IceFloeTracker, Random

julia> Random.seed!(123);

julia> bw_img = rand([0, 1], 5, 10)
5×10 Matrix{Int64}:
 1  0  1  0  0  0  0  0  0  1
 1  0  1  1  1  0  0  0  1  1
 1  1  0  1  1  0  1  0  0  1
 0  1  0  1  0  0  0  0  1  0
 1  0  0  0  0  1  0  1  0  1

julia> label_img = IceFloeTracker.label_components(bw_img, trues(3,3))
5×10 Matrix{Int64}:
 1  0  1  0  0  0  0  0  0  4
 1  0  1  1  1  0  0  0  4  4
 1  1  0  1  1  0  3  0  0  4
 0  1  0  1  0  0  0  0  4  0
 1  0  0  0  0  2  0  4  0  4

julia> properties = ["area", "perimeter"]
2-element Vector{String}:
 "area"
 "perimeter"

 julia> IceFloeTracker.regionprops_table(label_img, bw_img, properties = properties)
 4×2 DataFrame
  Row │ area   perimeter 
      │ Int32  Float64   
 ─────┼──────────────────
    1 │    13   11.6213
    2 │     1    0.0
    3 │     1    0.0
    4 │     7    4.62132
```
"""
function regionprops_table(
    label_img::Matrix{Int64},
    intensity_img::Union{Nothing,AbstractMatrix}=nothing;
    properties::Union{Vector{<:AbstractString},Tuple{String,Vararg{String}}}=(
        "centroid",
        "area",
        "major_axis_length",
        "minor_axis_length",
        "convex_area",
        "bbox",
        "perimeter",
    ),
    extra_properties::Union{Tuple{Function,Vararg{Function}},Nothing}=nothing,
)::DataFrame
    if !isnothing(extra_properties)
        @error "extra_properties not yet implemented in this wrapper; setting it to `nothing`"
        extra_properties = nothing
    end

    props = DataFrame(
        sk_measure.regionprops_table(
            label_img, intensity_img, properties; extra_properties=extra_properties
        ),
    )

    if "bbox" in properties
        bbox_cols = getbboxcolumns(props)
        fixzeroindexing!(props, bbox_cols[1:2])
        renamecols!(props, bbox_cols, ["min_row", "min_col", "max_row", "max_col"])
    end

    if "centroid" in properties
        centroid_cols = getcentroidcolumns(props)
        roundtoint!(props, centroid_cols)
        fixzeroindexing!(props, centroid_cols)
        renamecols!(props, centroid_cols, ["row_centroid", "col_centroid"])
    end
    return props
end

# Also adding regionprops as it might be more computationally efficient to get a property for each object on demand than generating the full regionprops_table at once

"""
    regionprops(label_img, ; properties, connectivity)

A wrapper of the `regionprops` function from the skimage python library.
    
See its full documentation at https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.regionprops.
    
# Arguments
- `label_img`: Image with the labeled objects of interest
- `intensity_img`: (Optional) Used for generating `extra_properties`, integer/float array from which (presumably) `label_img` was generated
- `extra_properties`: (Optional) not yet implemented. It will be set to `nothing`

See also [`regionprops_table`](@ref)

# Examples

```jldoctest; setup = :(using IceFloeTracker, Random)
julia> Random.seed!(123);

julia> bw_img = rand([0, 1], 5, 10)
5×10 Matrix{Int64}:
 1  0  1  0  0  0  0  0  0  1
 1  0  1  1  1  0  0  0  1  1
 1  1  0  1  1  0  1  0  0  1
 0  1  0  1  0  0  0  0  1  0
 1  0  0  0  0  1  0  1  0  1

 julia> label_img = Images.label_components(bw_img, trues(3,3))
 5×10 Matrix{Int64}:
  1  0  1  0  0  0  0  0  0  4
  1  0  1  1  1  0  0  0  4  4
  1  1  0  1  1  0  3  0  0  4
  0  1  0  1  0  0  0  0  4  0
  1  0  0  0  0  2  0  4  0  4

 julia> regions = regionprops(label_img, bw_img);

 julia> for region in regions
           println(region.area,"\t", region.perimeter)
        end
13      11.621320343559642
1       0.0
1       0.0
7       4.621320343559642
```
"""
function regionprops(
    label_img::Any,
    intensity_img::Any=nothing;
    extra_properties::Union{Tuple{Function,Vararg{Function}},Nothing}=nothing,
)
    if !isnothing(extra_properties)
        @error "<extra_properties> not yet implemented in this wrapper; setting it to <nothing>"
        extra_properties = nothing
    end

    return sk_measure.regionprops(
        label_img, intensity_img; extra_properties=extra_properties
    )
end

"""
    getcentroidcolumns(props::DataFrame)

Returns the column names of the centroid columns in `props`.
"""
function getcentroidcolumns(props::DataFrame)
    return filter(col -> occursin(r"^centroid-\d$", col), names(props))
end

"""
    getbboxcolumns(props::DataFrame)

Return the column names of the bounding box columns in `props`.
"""
function getbboxcolumns(props::DataFrame)
    return filter(col -> occursin(r"^bbox-\d$", col), names(props))
end

"""
    cropfloe(floesimg, props, i)

Crops the floe delimited by the bounding box data in `props` at index `i` from the floe image `floesimg`.
"""
function cropfloe(floesimg::BitMatrix, props::DataFrame, i::Int64)
    #= 
    Crop the floe using bounding box data in props.
    Note: Using a view of the cropped floe was considered but if there were multiple components in the cropped floe, the source array with the floes would be modified. =#
    prefloe = floesimg[props.min_row[i]:props.max_row[i], props.min_col[i]:props.max_col[i]]

    #= Check if more than one component is present in the cropped image.
    If so, keep only the largest component by removing all on pixels not in the largest component =#
    components = label_components(prefloe, trues(3, 3))

    if length(unique(components)) > 2
        mask = IceFloeTracker.bwareamaxfilt(components .> 0)
        prefloe[.!mask] .= 0
    end
    return prefloe
end

"""
    addfloearrays(props::DataFrame, floeimg::BitMatrix)

Add a column to `props` called `floearray` containing the cropped floe masks from `floeimg`.
"""
function addfloemasks!(props::DataFrame, floeimg::BitMatrix)
    props.mask = getfloemasks(props, floeimg)
    return nothing
end

"""
    getfloemasks(props::DataFrame, floeimg::BitMatrix)

Return a vector of cropped floe masks from `floeimg` using the bounding box data in `props`.
"""
function getfloemasks(props::DataFrame, floeimg::BitMatrix)
    return map(i -> cropfloe(floeimg, props, i), 1:nrow(props))
end

"""
    fixzeroindexing!(props::DataFrame, props_to_fix::Vector{T}) where T<:Union{Symbol,String}

Fix the zero-indexing of the `props_to_fix` columns in `props` by adding 1 to each element.
"""
function fixzeroindexing!(
    props::DataFrame, props_to_fix::Vector{T}
) where {T<:Union{Symbol,String}}
    props[:, props_to_fix] .+= 1
    return nothing
end

"""
    renamecols!(props::DataFrame, oldnames::Vector{T}, newnames::Vector{T}) where T<:Union{Symbol,String}

Rename the `oldnames` columns in `props` to `newnames`.
"""
function renamecols!(
    props::DataFrame, oldnames::Vector{T}, newnames::Vector{T}
) where {T<:Union{Symbol,String}}
    rename!(props, Dict(zip(oldnames, newnames)))
    return nothing
end

"""
    roundtoint!(props::DataFrame, colnames::Vector{T}) where T<:Union{Symbol,String}

Round the `colnames` columns in `props` to `Int`.
"""
function roundtoint!(props::DataFrame, colnames::Vector{T}) where {T<:Union{Symbol,String}}
    props[!, colnames] = (x -> convert.(Int, x))(round.(Int, props[!, colnames]))
    return nothing
end
