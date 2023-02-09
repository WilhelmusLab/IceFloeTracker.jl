"""
    extract_floe_features(bw::AbstractArray{Bool}, area_threshold::Tuple{Int64,Int64}, features::Union{Vector{Symbol},Vector{String}})::DataFrame

Extract features from the labeled image `bw` using the area thresholds in `area_threshold` and the vector of features `features`. Returns a DataFrame of the extracted features.

# Arguments
- `bw``: A labeled image.
- `area_threshold`: A tuple of the minimum and maximum (inclusive) areas of the floes to extract.
- `features`: A vector of the features to extract.

# Example
```jldoctest; setup = :(using IceFloeTracker, Random)
julia> Random.seed!(123)
TaskLocalRNG()

julia> bw_img = rand(Bool, 5, 50)
5×50 Matrix{Bool}:
 0  0  1  0  1  0  1  0  1  0  1  0  0  0  0  0  0  1  1  0  1  0  1  1  1  1  1  0  0  0  0  0  0  0  1  1  1  0  0  1  1  1  0  0  1  1  1  1  0  1
 1  0  0  0  1  1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  1  0  0  1  0  0  0  0  1  1  1  1  1  0  1  0  1  1  0  1  1  1  1  0  0  1  1  1  1
 0  1  1  1  1  1  0  0  1  0  1  1  0  1  1  1  1  1  1  1  1  0  0  1  0  1  1  1  0  0  0  0  1  0  0  0  0  0  1  0  0  0  1  0  1  0  0  0  0  0
 0  1  0  0  1  1  1  0  0  0  1  0  1  0  1  0  1  0  0  1  1  1  1  0  0  1  0  0  0  0  1  0  0  0  1  1  1  0  0  1  1  1  0  1  0  0  0  0  1  0
 0  1  1  1  1  0  0  0  1  0  1  0  0  1  0  1  0  0  1  1  0  1  1  0  1  1  0  1  1  0  0  0  1  0  0  0  0  1  0  0  0  1  1  0  1  0  1  1  1  0

julia> properties = ["centroid", "area", "major_axis_length", "minor_axis_length", "convex_area", "bbox"]
6-element Vector{String}:
 "centroid"
 "area"
 "major_axis_length"
 "minor_axis_length"
 "convex_area"
 "bbox"

julia> area_threshold = (1, 5)
(1, 5)

julia> feats = IceFloeTracker.extract_floe_features(bw_img, area_threshold, properties)
8×10 DataFrame
 Row │ area   bbox-0  bbox-1  bbox-2  bbox-3  centroid-0  centroid-1  convex_area  major_axis_length  minor_axis_length 
     │ Int32  Int32   Int32   Int32   Int32   Float64     Float64     Int32        Float64            Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │     1       1       3       2       4        0.0         2.0             1            0.0                0.0
   2 │     1       5       9       6      10        4.0         8.0             1            0.0                0.0
   3 │     2       1      18       2      20        0.0        17.5             2            2.0                0.0
   4 │     2       5      28       6      30        4.0        27.5             2            2.0                0.0
   5 │     1       4      31       5      32        3.0        30.0             1            0.0                0.0
   6 │     1       5      33       6      34        4.0        32.0             1            0.0                0.0
   7 │     4       4      35       6      39        3.25       35.5             5            4.68021            1.04674
   8 │     4       4      47       6      50        3.75       47.25            5            3.4641             1.41421
 
```
"""
function extract_floe_features(bw::T; area_threshold::Tuple{Int64,Int64}=(300, 90000), features::Union{Vector{Symbol},Vector{<:AbstractString}})::DataFrame where {T<:AbstractArray{Bool}}
    # assert the first area threshold is less than the second
    check_2_tuple(area_threshold)

    props = regionprops_table(label_components(bw, trues(3, 3)), properties=features)

    # filter by area using the area threshold
    return props[area_threshold[1].<=props.area.<=area_threshold[2], :]

    # TODO: floe post-processing
end

function extract_floe_features(; input::String, output::String, area_threshold::String, features::String)::Vector{DataFrame}
    # parse the area threshold
    area_threshold = parse_2tuple(area_threshold)

    # parse the features
    features = split(features) # need to check this works

    # assert the first area threshold is less than the second
    check_2_tuple(area_threshold)

    # load segmented images in input directory
    segmented_floes = [
        BitMatrix(load(joinpath(input, f))) for
        f in readdir(input)
    ]

    props = [
        IceFloeTracker.extract_floe_features(bw; area_threshold=area_threshold, features=features) for
        bw in segmented_floes]

    # serialize the props vector to the output directory using JLD2
    serialize(joinpath(output, "floe_library.dat"), props) # need job id for file name?
    return props
end