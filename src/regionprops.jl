
"""
    regionprops_table(label_img, intensity_img; properties, connectivity, extra_properties)

A wrapper of the `regionprops_table` function from the skimage python library.
    
See its full documentation at https://scikit-image.org/docs/stable/api/skimage.measure.html#regionprops-table.
    
# Arguments
- `label_img`: Image with the labeled objects of interest
- `intensity_img`: (Optional) Used for generating `extra_properties`, integer/float array from which (presumably) `label_img` was generated 
- `properties`: List (`Vector` or `Tuple`) of properties to be generated for each connected component in `label_img`
- `extra_properties`: (Optional) not yet implemented. It will be set to `nothing`

See also [`regionprops`](@ref)

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

julia> properties = ["area", "perimeter"]
2-element Vector{String}:
 "area"
 "perimeter"

 julia> regionprops_table(label_img, bw_img, properties = properties, dataframe=true)
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
    label_img::Any,
    intensity_img::Any=nothing;
    properties::Union{Vector{String},Tuple{String,Vararg{String}}}=("centroid", "area", "major_axis_length", "minor_axis_length", "convex_area", "bbox"),
    extra_properties::Union{Tuple{Function,Vararg{Function}},Nothing}=nothing,
)::DataFrame
    if !isnothing(extra_properties)
        @error "extra_properties not yet implemented in this wrapper; setting it to `nothing`"
        extra_properties = nothing
    end

    props = sk_measure.regionprops_table(
        label_img, intensity_img, properties; extra_properties=extra_properties
    ) |> DataFrame

    # Add one to bbox-* cols to account for 0-based indexing
    if "bbox" in properties
        props[:,["bbox-0", "bbox-1", "bbox-2", "bbox-3"]] .+= 1
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
