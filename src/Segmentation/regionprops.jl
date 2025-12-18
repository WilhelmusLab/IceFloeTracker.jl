
import ..skimage: sk_measure
import DataFrames: rename!, DataFrame, nrow, select!
import ..Geospatial: latlon
import Images: 
    component_boxes, 
    component_centroids, 
    component_lengths, 
    erode,
    Fill,
    labels_map, 
    padarray,
    SegmentedImage,
    strel_diamond,
    strel_box

import DSP: conv

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

```jldoctest; setup = :(using IceFloeTracker, Random, Images)
julia> using IceFloeTracker, Random, Images

julia> Random.seed!(123);

julia> bw_img = rand([0, 1], 5, 10)
5×10 Matrix{Int64}:
 1  0  1  0  0  0  0  0  0  1
 1  0  1  1  1  0  0  0  1  1
 1  1  0  1  1  0  1  0  0  1
 0  1  0  1  0  0  0  0  1  0
 1  0  0  0  0  1  0  1  0  1

julia> label_img = label_components(bw_img, trues(3,3))
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

 julia> regionprops_table(label_img, bw_img, properties = properties)
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
# function regionprops_table(
#     label_img::Union{Matrix{Int64},SegmentedImage},
#     intensity_img::Union{Nothing,AbstractMatrix}=nothing;
#     properties::Union{Vector{<:AbstractString},Tuple{String,Vararg{String}}}=(
#         "label",
#         "centroid",
#         "area",
#         "major_axis_length",
#         "minor_axis_length",
#         "convex_area",
#         "bbox",
#         "perimeter",
#         "orientation",
#     ),
#     extra_properties::Union{Tuple{Function,Vararg{Function}},Nothing}=nothing,
# )::DataFrame
#     if label_img isa SegmentedImage
#         label_img = labels_map(label_img)
#     end

#     if !isnothing(extra_properties)
#         @error "extra_properties not yet implemented in this wrapper; setting it to `nothing`"
#         extra_properties = nothing
#     end

#     props = DataFrame(
#         sk_measure.regionprops_table(
#             label_img, intensity_img, properties; extra_properties=extra_properties
#         ),
#     )

#     if "bbox" in properties
#         bbox_cols = getbboxcolumns(props)
#         fixzeroindexing!(props, bbox_cols[1:2])
#         renamecols!(props, bbox_cols, ["min_row", "min_col", "max_row", "max_col"])
#     end

#     if "centroid" in properties
#         centroid_cols = getcentroidcolumns(props)
#         roundtoint!(props, centroid_cols)
#         fixzeroindexing!(props, centroid_cols)
#         renamecols!(props, centroid_cols, ["row_centroid", "col_centroid"])
#     end
#     return props
# end

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

"""
    addlatlon(pairedfloesdf::DataFrame, refimage::AbstractString)

Add columns `latitude`, `longitude`, and pixel coordinates `x`, `y` to `pairedfloesdf`.

# Arguments
- `pairedfloesdf`: dataframe containing floe tracking data.
- `refimage`: path to reference image.
"""
function addlatlon!(pairedfloesdf::DataFrame, refimage::AbstractString)
    latlondata = latlon(refimage)
    colstodrop = [:row_centroid, :col_centroid, :min_row, :min_col, :max_row, :max_col]
    converttounits!(pairedfloesdf, latlondata, colstodrop)
    return nothing
end

## LatLon functions originally from IFTPipeline.jl

"""
    convertcentroid!(propdf, latlondata, colstodrop)

Convert the centroid coordinates from row and column to latitude and longitude dropping unwanted columns specified in `colstodrop` for the output data structure. Addionally, add columns `x` and `y` with the pixel coordinates of the centroid.
"""
function convertcentroid!(propdf, latlondata, colstodrop)
    latitude, longitude = [
        [
            latlondata[c][Int(round(x)), Int(round(y))] for
            (x, y) in zip(propdf.row_centroid, propdf.col_centroid)
        ] for c in [:latitude, :longitude]
    ]

    x, y = [
        [latlondata[c][Int(round(z))] for z in V] for
        (c, V) in zip([:Y, :X], [propdf.row_centroid, propdf.col_centroid])
    ]

    propdf.latitude = latitude
    propdf.longitude = longitude
    propdf.x = x
    propdf.y = y
    select!(propdf, Not(colstodrop))
    return nothing
end

"""
    converttounits!(propdf, latlondata, colstodrop)

Convert the floe properties from pixels to kilometers and square kilometers where appropiate. Also drop the columns specified in `colstodrop`.
"""
function converttounits!(propdf, latlondata, colstodrop)
    if nrow(propdf) == 0
        select!(propdf, Not(colstodrop))
        insertcols!(
            propdf,
            :latitude => Float64,
            :longitude => Float64,
            :x => Float64,
            :y => Float64,
        )
        return nothing
    end
    convertcentroid!(propdf, latlondata, colstodrop)
    x = latlondata[:X]
    dx = abs(x[2] - x[1])
    convertarea(area) = area * dx^2 / 1e6
    convertlength(length) = length * dx / 1e3
    propdf.area .= convertarea(propdf.area)
    propdf.convex_area .= convertarea(propdf.convex_area)
    propdf.minor_axis_length .= convertlength(propdf.minor_axis_length)
    propdf.major_axis_length .= convertlength(propdf.major_axis_length)
    propdf.perimeter .= convertlength(propdf.perimeter)
    return nothing
end

"""component_floes(indexmap; minimum_area=1)

Given a labeled array, produce a dictionary where each entry is a binary array cropped to the
location of the label.

Examples:
```jldoctest; setup = :(using IceFloeTracker)

julia> A = zeros(Int, 8, 6); A[2:6, 1:2] .= 1; A[3:7, 4:5] .= 2; A[4:6, 3:6] .= 2;
julia> A
8×6 Matrix{Int64}:
 0  0  0  0  0  0
 1  1  0  0  0  0
 1  1  0  2  2  0
 1  1  2  2  2  2
 1  1  2  2  2  2
 1  1  2  2  2  2
 0  0  0  2  2  0
 0  0  0  0  0  0

julia> masks = component_floes(A)
julia> Dict{Any, Any} with 3 entries:
  0 => Bool[1 1 … 1 1; 0 0 … 1 1; … ; 1 1 … 0 1; 1 1 … 1 1]
  2 => Bool[0 1 1 0; 1 1 1 1; … ; 1 1 1 1; 0 1 1 0]
  1 => Bool[1 1; 1 1; … ; 1 1; 1 1]
julia> masks[2]
5×4 BitMatrix:
 0  1  1  0
 1  1  1  1
 1  1  1  1
 1  1  1  1
 0  1  1  0
julia> masks[1]
5×2 BitMatrix:
 1  1
 1  1
 1  1
 1  1
 1  1
```
"""
function component_floes(indexmap; minimum_area=1)
    labels = unique(indexmap)
    boxes = component_boxes(indexmap)
    mn, mx = extrema(indexmap)
    if !(mn == 0 || mn == 1)
        throw(ArgumentError("The input labeled array should contain background label `0` as the minimum value"))
    end
    areas = component_lengths(indexmap)
    floe_masks = Dict()
    for i in labels
        if areas[i] > minimum_area
            floe_masks[i] = indexmap[boxes[i]] .== i
        end
    end
    return floe_masks
end

"""component_perimeter(indexmap; minimum_area=1, algorithm="benkrid_crookes")

Estimate the perimeter of the labeled regions in `indexmap` using the specified algorithm.
Algorithm options = "benkrid_crookes" (only option currently, will add crofton in future release)


"""
function component_perimeters(
    indexmap;
    minimum_area=1, 
    algorithm="benkrid_crookes",
    connectivity=4)
    masks = component_floes(indexmap)
    perims = Dict()
    algorithm ∉ ["benkrid_crookes"] && begin
        print("Unsupported Algorithm, defaulting to Benkrid-Crookes")
        algorithm = "benkrid_crookes"
    end

    algorithm_functions = Dict("benkrid_crookes" => benkrid_crookes)

    for label in keys(masks)
        n, m = size(masks[label])
        n * m == 1 && continue
        label == 0 && (masks[label] = 0; continue)
        
        # Shape needs to have a border of zeros for erode to work here    
        mpad = padarray(masks[label], Fill(0, (1, 1)))
        connectivity == 4 ? (strel = strel_diamond((3,3))) : (strel = strel_box((3,3)))
        epad = mpad .- erode(mpad, strel)
        e = epad[1:n, 1:m]
        perims[label] = algorithm_functions[algorithm](e)
    end
    return perims
end

function benkrid_crookes(edge_array)
    type_vals = zeros(33)
    type_vals[[5, 7, 15, 17, 25, 27]] .= 1
    type_vals[[21, 33]] .= sqrt(2)
    type_vals[[13, 23]] .= (1 + sqrt(2)) / 2

    conv_arr = [10 2 10; 2 1 2; 10 2 10]
    results = conv(edge_array, conv_arr; algorithm=:direct)

    val_counts = Dict()
    for val in vec(results)
        val ∉ keys(val_counts) && (val_counts[val] = 0)
        val_counts[val] += 1
    end
    perim = 0
    for val in keys(val_counts)
        val == 0 && continue
        val > 33 && continue
        perim += type_vals[val] * val_counts[val]
    end
    
    return perim
end


# TODO: Make the new function work with the old docstring and example intact
"""
    regionprops_table(label_img, intensity_img; properties, connectivity, extra_properties)

A wrapper of the `regionprops_table` function from the skimage python library.
    
See its full documentation at https://scikit-image.org/docs/stable/api/skimage.measure.html#regionprops-table.
    
# Arguments
- `label_img`: Image with the labeled objects of interest. May be an integer array or a SegmentedImage.
- `intensity_img`: (Optional) Used for generating `extra_properties`, such as a color image to use for calculating mean color in segments.
- `properties`: List (`Vector` or `Tuple`) of properties to be generated for each connected component in `label_img`
- `extra_properties`: (Optional) not yet implemented. It will be set to `nothing`

# Notes
- Zero indexing has been corrected for the `bbox` and `centroid` properties
- `bbox` data (`max_col` and `max_row`) are inclusive
- `centroid` data are rounded to the nearest integer

See also [`regionprops`](@ref)

# Examples

```jldoctest; setup = :(using IceFloeTracker, Random, Images)
julia> using IceFloeTracker, Random, Images

julia> Random.seed!(123);

julia> bw_img = rand([0, 1], 5, 10)
5×10 Matrix{Int64}:
 1  0  1  0  0  0  0  0  0  1
 1  0  1  1  1  0  0  0  1  1
 1  1  0  1  1  0  1  0  0  1
 0  1  0  1  0  0  0  0  1  0
 1  0  0  0  0  1  0  1  0  1

julia> label_img = label_components(bw_img, trues(3,3))
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

 julia> regionprops_table(label_img, bw_img, properties = properties)
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
    label_img::Union{Matrix{Int64},SegmentedImage},
    intensity_img::Union{Nothing,AbstractMatrix}=nothing;
    properties::Union{Vector{<:AbstractString}, Vector{<:Symbol}}=[
        :label,
        :centroid,
        :area,
        :major_axis_length,
        :minor_axis_length,
        :convex_area,
        :bbox,
        :perimeter,
        :orientation,
    ],
    extra_properties::Union{Tuple{Function,Vararg{Function}},Nothing}=nothing,
    minimum_area=1,
)::DataFrame

    typeof(label_img) == SegmentedImage ? (labels = labels_map(label_img)) : labels = deepcopy(label_img)
    eltype(properties) == String && (properties = [Symbol(a) for a in properties])

    data = Dict{Symbol,Any}()
    
    # Begin by extracting the set of labels that meet the minimum area criterion
    # We also get a sorted list of image labels, so all the dictionary entries can 
    # be placed in the same order.
    areas = segment_pixel_count(label_img)
    img_labels = segment_labels(label_img)
    img_labels = img_labels[img_labels .!= 0]
    sort!(img_labels)
    img_labels = img_labels[[areas[s] > minimum_area for s in img_labels]]

    :label ∈ properties && push!(data, :label => img_labels)

    # re-use the computed areas
    :area ∈ properties && begin
        push!(data, :area => map(s -> areas[s], img_labels))
    end

    :mask ∈ properties && begin
        floe_masks = component_floes(labels_map(label_img); minimum_area=minimum_area)
        push!(data, :mask => map(s -> floe_masks[s], img_labels))
    end

    :perimeter ∈ properties && begin
        # future option: allow component perimeters to take the floe masks as an argument to save compute time
        floe_perims = component_perimeters(labels_map(label_img); minimum_area=minimum_area)
        push!(data, :perimeter => map(s -> floe_perims[s], img_labels)) 
    end


    # TODO measurements
    # get centroid and convert to row_centroid, col_centroid
    # get major/minor axis
    # get orientation
    # get perimeter
    # get convex_area
    # bounding boxes

    return DataFrame(data)[:, properties]
end