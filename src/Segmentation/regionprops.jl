import DataFrames: rename!, DataFrame, nrow, select!
import ..Geospatial: latlon
import Images: 
    component_boxes, 
    component_centroids, 
    component_lengths, 
    component_indices,
    convexhull,
    erode,
    Fill,
    labels_map, 
    padarray,
    SegmentedImage,
    strel_diamond,
    strel_box

import DSP: conv

abstract type PerimeterEstimationAlgorithm <: Function end
abstract type ConvexAreaEstimationAlgorithm <: Function end


# TODO: Determine if this function is needed, since rename! is a standard function for DataFrames.jl.
"""
    renamecols!(props::DataFrame, oldnames::Vector{T}, newnames::Vector{T}) where T<:Union{Symbol,String}

Rename the `oldnames` columns in `props` to `newnames`.
"""
function renamecols!(
    props::DataFrame,
    oldnames::Vector{Union{Symbol,String}},
    newnames::Vector{Union{Symbol,String}}
    )
    rename!(props, Dict(zip(oldnames, newnames)))
    return nothing
end

# TODO: Check if we can drop this function. It's actually useful to keep the decimals in the row/col centroids.
"""
    roundtoint!(props::DataFrame, colnames::Vector{T}) where T<:Union{Symbol,String}

Round the `colnames` columns in `props` to `Int`.
"""
function roundtoint!(
    props::DataFrame,
    colnames::Vector{Union{Symbol,String}}
    )
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
function addlatlon!(
    pairedfloesdf::DataFrame,
    refimage::AbstractString
    )
    latlondata = latlon(refimage)
    colstodrop = [:row_centroid, :col_centroid, :min_row, :min_col, :max_row, :max_col]
    converttounits!(pairedfloesdf, latlondata, colstodrop)
    return nothing
end

## LatLon functions originally from IFTPipeline.jl
# TODO: Add example with reference geotiff image.
"""
    convertcentroid!(propdf, latlondata, colstodrop)

Convert the centroid coordinates from row and column to latitude and longitude dropping unwanted
columns specified in `colstodrop` for the output data structure. Addionally, add columns `x` and `y`
with the pixel coordinates of the centroid.
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
    mn = minimum(indexmap)
    if !(mn == 0 || mn == 1)
        throw(ArgumentError("The input labeled array should contain background label `0` as the minimum value"))
    end
    areas = component_lengths(indexmap)
    floe_masks = Dict(i => indexmap[boxes[i]] .== i for i in labels if areas[i] > minimum_area)
    return floe_masks
end

"""component_perimeters(indexmap; minimum_area=1, algorithm="benkrid_crookes")

Estimate the perimeter of the labeled regions in `indexmap` using the specified algorithm.
Algorithm options = "benkrid_crookes" (only option currently, will add crofton in future release)
Defaults to using connectivity 4.
"""
function component_perimeters(
    indexmap;
    algorithm::PerimeterEstimationAlgorithm=BenkridCrookes(),
    )
    masks = component_floes(indexmap)
    perims = Dict(
        label => (label == 0 ? 0 : algorithm(masks[label])) for
        label in keys(masks) if prod(size(masks[label])) > 1
    )
    return perims
end

# TODO: Test and implement correction factor (multiply B-K perimeter by 0.95 if larger than some factor.)
# TODO: Implement the crofton perimeter algorithm
"""BenkridCrookes(connectivity=4)
   
Functor producing a BenkridCrookes PerimeterEstimationAlgorithm. The connectivity
used for the erosion is the only parameter. The algorithm uses strel_diamond((3,3)) for 4-connectivity and
strel_box((3,3)) for 8-connectivity. The resulting function operates on a binary array, which 
is assumed to contain a single object.

# Examples
```
julia> A = [0 1 1; 1 1 1; 1 1 1];

julia> BenkridCrookes(connectivity=4)(A)
7.414213562373095

julia> BenkridCrookes(connectivity=8)(A)
7.0
```
"""
@kwdef struct BenkridCrookes <: PerimeterEstimationAlgorithm
    connectivity = 4
end

function (f::BenkridCrookes)(shape_array)
    f.connectivity == 4 ? (strel = strel_diamond((3,3))) : (strel = strel_box((3,3)))
    # Get border using the strel
    # Shape needs to have a border of zeros for erode to work here    
    n, m = size(shape_array)
    mpad = padarray(shape_array, Fill(0, (1, 1)))    
    epad = mpad .- erode(mpad, strel)
    e = epad[1:n, 1:m]

    # Set up lookup table for computing perimeter
    type_vals = zeros(33)
    type_vals[[5, 7, 15, 17, 25, 27]] .= 1
    type_vals[[21, 33]] .= sqrt(2)
    type_vals[[13, 23]] .= (1 + sqrt(2)) / 2

    # Convolution array for classifying boundary pixel type
    conv_arr = [10 2 10; 2 1 2; 10 2 10]

    results = conv(e, conv_arr; algorithm=:direct)

    # Count instances of boundary types and multiply to get the perimeter
    val_counts = Dict()
    for val in vec(results)
        val_counts[val] = get(val_counts, val, 0) + 1
    end
    perim = sum(
        type_vals[val] * count for (val, count) in pairs(val_counts) if val > 0 && val <= 33
    )
    
    return perim
end


"""component_convex_area(A; algorithm=PixelConvexArea()")

Compute the convex area of labeled regions. Two methods available: "PixelConvexArea" and "PolygonConvexArea".
The polygon method uses Green's theorem to find the area of a polygon through its line integral, 
while the pixel method uses a point-in-pixel calculation to determine if pixels are inside the
convex hull. In general the polygon area will be smaller than the pixel area.
"""
function component_convex_areas(A;
    algorithm::ConvexAreaEstimationAlgorithm=PixelConvexArea()
    ) 
    mn = minimum(A)

    (mn != 0 && mn != 1) && throw(
        ArgumentError(
            "The input labeled array should contain background label `0` as the minimum value",
        ),
    )

    return algorithm(A)
end

"""PolygonConvexArea(minimum_area=4)

Estimate the convex area by integrating the area of the convex hull polygon.
Uses the convexhull function from ImageMorphology, which raises an error if the
area of the segment is less than or equal to 3. In general, the error should be smaller
for larger shapes.
"""
@kwdef struct PolygonConvexArea <: ConvexAreaEstimationAlgorithm
    minimum_area = 4
end

function (f::PolygonConvexArea)(A)
    mx = maximum(A)
    areas = component_lengths(A)

    convex_areas = zeros(Float64, 0:mx)
    for i in unique(A)
        # treat convex area background and too-small objects as undefined
        (i == 0) || (areas[i] < f.minimum_area) && begin
            convex_areas[i] = NaN
            continue
        end

        chull = convexhull(A .== i)
        N = length(chull)

        ca = 0
        for j in 1:N
            x0, y0 = Tuple(chull[j])
            x1, y1 = Tuple(chull[(j % N) + 1])
            ca += x0*y1 - y0*x1
        end
        ca *= 0.5
        convex_areas[i] = ca
    end
    return convex_areas
end


"""PixelConvexArea(minimum_area)
   PixelConvexArea(A)

Estimate the convex area by summing the number of pixels with indices falling within the
polygon formed by the convex hull. Uses the convexhull function from ImageMorphology, which raises an error if the
area of the segment is less than or equal to 3. In general, the error should be smaller
for larger shapes.
"""
@kwdef struct PixelConvexArea <: ConvexAreaEstimationAlgorithm
    minimum_area = 4
end

function (f::PixelConvexArea)(A)
    mx = maximum(A)
    convex_areas = zeros(Float64, 0:mx)
    areas = component_lengths(A)
    bboxes = component_boxes(A)
    labels = unique(A)
    for i in labels
        
        # treat convex area background and too-small objects as undefined
        (i == 0) || (areas[i] < f.minimum_area) && begin
            convex_areas[i] = NaN
            continue
        end
            
        chull = convexhull(A .== i)
        N = length(chull)
        x = getindex.(bboxes[i], 1)
        y = getindex.(bboxes[i], 2)
        
        for idx in eachindex(x)
            xi, yi = x[idx], y[idx]
            A[xi, yi] .== i && (convex_areas[i] += 1, continue)
            checkvals = zeros(N)
            for j in 1:N
                x0, y0 = Tuple(chull[j])
                x1, y1 = Tuple(chull[(j % N) + 1])
                checkvals[j] = (yi - y0)*(x1 - x0) - (xi - x0)*(y1 - y0)
            end
            all(checkvals .>= 0) && (convex_areas[i] += 1) 
            end
        end    
    return convex_areas
end

"""
    regionprops_table(label_img, intensity_img; properties, connectivity, extra_properties)

Compute measures of labeled regions in label_img and return as a DataFrame. Optionally, include an
extra image or array associated with the labels.
        
# Arguments
- `label_img`: Image with the labeled objects of interest. May be an integer array or a SegmentedImage.
- `intensity_img`: (Optional) Used for generating `extra_properties`, such as a color image to use for calculating mean color in segments.
- `properties`: List (`Vector` or `Tuple`) of properties to be generated for each connected component in `label_img`
- `extra_properties`: (Optional) not yet implemented. It will be set to `nothing`
- `minimum_area`: Smallest region to calculate measures on
- `perimeter_algorithm`: PerimeterEstimationAlgorithm. Currently only available is BenkridCrookes.
- `convex_area_algorithm`: ConvexAreaEstimationAlgorithm. Options are PolygonConvexArea or PixelConvexArea, default is pixel.

# Notes
- `bbox` data (`max_col` and `max_row`) are inclusive

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
    perimeter_algorithm=BenkridCrookes(),
    convex_area_algorithm=PixelConvexArea()
    )::DataFrame

    data = regionprops(label_img, intensity_img;
        properties=properties,
        extra_properties=extra_properties,
        minimum_area=minimum_area,
        perimeter_algorithm=perimeter_algorithm,
        convex_area_algorithm=convex_area_algorithm)
    return DataFrame(data)
end

# TODO: Updated regionprops_table and regionprops to enable customizing the algorithms.

"""regionprops(label_img, intensity_img; properties, extra_properties, minimum_area)

Core function returning a dictionary with an entry for each returned property.

See also [`regionprops_table`](@ref)

# Arguments
- `label_img`: Image with the labeled objects of interest. May be an integer array or a SegmentedImage.
- `intensity_img`: (Optional) Used for generating `extra_properties`, such as a color image to use for calculating mean color in segments.
- `properties`: List (`Vector` or `Tuple`) of properties to be generated for each connected component in `label_img`
- `extra_properties`: (Optional) not yet implemented. It will be set to `nothing`
- `minimum_area`: Smallest region to calculate measures on
- `perimeter_algorithm`: PerimeterEstimationAlgorithm. Currently only available is BenkridCrookes.
- `convex_area_algorithm`: ConvexAreaEstimationAlgorithm. Options are PolygonConvexArea or PixelConvexArea, default is pixel.


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

 julia> measures = regionprops(label_img, bw_img);

 julia> for s in unique(label_img)
           println(measures.area[s],"\t", measures.perimeter[s])
        end
13      11.621320343559642
1       0.0
1       0.0
7       4.621320343559642
```
"""
function regionprops(
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
    perimeter_algorithm=BenkridCrookes(),
    convex_area_algorithm=PixelConvexArea()
)

    isa(label_img, SegmentedImage) ? (labels = labels_map(label_img)) : labels = label_img
    eltype(properties) == String && (properties = [Symbol(a) for a in properties])

    data = Dict{Symbol,Any}()
    
    # Begin by extracting the set of labels that meet the minimum area criterion
    # We also get a sorted list of image labels, so all the dictionary entries can 
    # be placed in the same order.
    areas = component_lengths(labels)
    img_labels = unique(labels)
    img_labels = img_labels[img_labels .!= 0]
    sort!(img_labels)
    img_labels = img_labels[[areas[s] > minimum_area for s in img_labels]]

    :label ∈ properties && push!(data, :label => img_labels)

    # re-use the computed areas
    :area ∈ properties && begin
        push!(data, :area => map(s -> areas[s], img_labels))
    end

    # Centroid, major/minor axis, and orientation all come from the image moments.
    compute_moments = any([:centroid ∈ properties,
                          :major_axis_length ∈ properties,
                          :minor_axis_length ∈ properties,
                          :orientation ∈ properties])
    compute_moments && begin
        data_moments = _component_moment_measures(labels, img_labels)
        :centroid ∈ properties && begin
            push!(data, :row_centroid => data_moments[:row_centroid])
            push!(data, :col_centroid => data_moments[:col_centroid])
        end 
        :major_axis_length ∈ properties && push!(data, :major_axis_length => data_moments[:major_axis_length])
        :minor_axis_length ∈ properties && push!(data, :minor_axis_length => data_moments[:minor_axis_length])
        :orientation ∈ properties  && push!(data, :orientation => data_moments[:orientation])
    end

    :bbox ∈ properties && begin
        bboxes_init =  component_boxes(labels)
        bboxes = stack(_get_bounds.(bboxes_init[s] for s in img_labels))

        push!(data, :min_row => bboxes[1,:]) 
        push!(data, :max_row => bboxes[2,:]) 
        push!(data, :min_col => bboxes[3,:]) 
        push!(data, :max_col => bboxes[4,:]) 
    end 

    :perimeter ∈ properties && begin
        floe_perims = component_perimeters(labels; algorithm=perimeter_algorithm)
        push!(data, :perimeter => map(s -> floe_perims[s], img_labels)) 
    end

    :convex_area ∈ properties && begin
        convex_areas = component_convex_areas(labels; algorithm=convex_area_algorithm)
        push!(data, :convex_area => map(s -> convex_areas[s], img_labels))
    end


    # psi-s needs masks, so this can get called first
    :mask ∈ properties && begin
        floe_masks = component_floes(labels; minimum_area=minimum_area)
        push!(data, :mask => map(s -> floe_masks[s], img_labels))
    end

    # TODO add psi-s curve generation as an option
    updated_properties = Symbol[]
    for p in properties
        if p == :bbox
            append!(updated_properties, [:min_row, :max_row, :min_col, :max_col])
        elseif p == :centroid
            append!(updated_properties, [:row_centroid, :col_centroid])
        else
            append!(updated_properties, [p])
        end
    end
    return Dict(prop => data[prop] for prop in updated_properties)
end

"""_get_bounds(box)

Helper function to extract the min row, max row, min col, and max col
from a vector of CartesianIndex ranges
"""
function _get_bounds(box)
    upper_left = first(collect(box))
    lower_right = last(collect(box))
    minrow = getindex(upper_left, 1)
    mincol = getindex(upper_left, 2)
    maxrow = getindex(lower_right, 1)
    maxcol = getindex(lower_right, 2)
    return minrow, maxrow, mincol, maxcol
end

"""_central_moment(x, y, xc, yc, p, q)
Compute the central moments based on the index vectors
x and y, the centroid xc, yc, and the exponents p and q.
"""
function _central_moment(x, y, xc, yc, p, q)
    mu = sum((x .- xc).^p .* (y .- yc).^q)
    return mu
end

"""_component_moment_measures(labels, label_list)

Compute standard image measures which rely on calculation of image moments. Returns
a dictionary

# Arguments
- labels: Matrix with integer-labeled regions
- label_list: List of labels to compute measures for

# Outputs
- Dictionary with entries "row_centroid", "col_centroid", "major_axis_length", "minor_axis_length", 
and "orientation" where each entry is a vector ordered by label_list.
"""
function _component_moment_measures(labels, label_list)
    data = Dict()
    centroids = component_centroids(labels)
    areas = component_lengths(labels)
    indices = component_indices(CartesianIndex, labels)

    row_centroid = first.(centroids)
    col_centroid = last.(centroids)
    push!(data, :row_centroid => map(s -> row_centroid[s], label_list))
    push!(data, :col_centroid => map(s -> col_centroid[s], label_list))

    moment_measures = []
    for s in label_list
        X = getindex.(indices[s], 1)
        Y = getindex.(indices[s], 2);
        xc = row_centroid[s]
        yc = col_centroid[s]
        
        μ11 = _central_moment(X, Y, xc, yc, 1, 1)
        μ20 = _central_moment(X, Y, xc, yc, 2, 0)
        μ02 = _central_moment(X, Y, xc, yc, 0, 2)
        θ = 0.5 * atan(2*μ11, μ20 - μ02)
        
        λ0 = (μ20 + μ02 + sqrt((μ20 - μ02)^2 + 4*μ11^2))/2
        λ1 = (μ20 + μ02 - sqrt((μ20 - μ02)^2 + 4*μ11^2))/2
        
        ra = 4*(λ0/areas[s])^0.5 
        rb = 4*(λ1/areas[s])^0.5
        append!(moment_measures, [[ra, rb, θ]])
    end
    moment_measures = stack(moment_measures)
    push!(data, :major_axis_length => moment_measures[1,:])
    push!(data, :minor_axis_length => moment_measures[2,:])
    push!(data, :orientation => moment_measures[3,:])

    return data
end
