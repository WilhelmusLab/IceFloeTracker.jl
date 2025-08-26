"""
    kmeans_segmentation(gray_image, ice_labels;)

Apply k-means segmentation to a gray image to isolate a cluster group representing sea ice. Returns a binary image with ice segmented from background.

# Arguments
- `gray_image`: output image from `ice-water-discrimination.jl` or gray ice floe leads image in `segmentation_f.jl`
- `ice_labels`: vector if pixel coordinates output from `find_ice_labels.jl`

"""
function kmeans_segmentation(
    gray_image::Matrix{Gray{Float64}},
    ice_labels::Union{Vector{Int64},BitMatrix,AbstractArray{<:Gray}},
    k::Int64=4,
    maxiter::Int64=50,
)::BitMatrix
    segmented = kmeans_segmentation(gray_image; k=k, maxiter=maxiter)

    ## Isolate ice floes and contrast from background
    segmented_ice = get_segmented_ice(segmented, ice_labels)
    return segmented_ice
end

function kmeans_segmentation(
    gray_image::Matrix{Gray{Float64}}; k::Int64=4, maxiter::Int64=50, random_seed::Int64=45
)
    Random.seed!(random_seed)

    ## NOTE(tjd): this clusters into k classes and solves iteratively with a max of maxiter iterations
    feature_classes = Clustering.kmeans(
        vec(gray_image), k; maxiter=maxiter, display=:none, init=:kmpp
    )

    class_assignments = assignments(feature_classes)

    segmented = reshape(class_assignments, size(gray_image))

    return segmented
end

function get_segmented_ice(
    segmented::Matrix{Int64}, ice_labels::Union{Vector{Int64},BitMatrix}
)
    ## Same principle as the get_nlabels function
    ## Has the weakness that only one segment can ever be chosen.
    isempty(ice_labels) && return falses(size(segmented))
    return segmented .== StatsBase.mode(segmented[ice_labels])
end

function get_segmented_ice(segmented::Matrix{Int64}, ice_labels::AbstractArray{<:Gray})
    boolean_map = ice_labels .|> Bool
    !(any(boolean_map)) && return falses(size(segmented))
    return segmented .== StatsBase.mode(segmented[boolean_map])
end

"""
    segmented_ice_cloudmasking(gray_image, cloudmask, ice_labels;)

Apply cloudmask to a bitmatrix of segmented ice after kmeans clustering. Returns a bitmatrix with open water/clouds = 0, ice = 1).

# Arguments

- `gray_image`: output image from `ice-water-discrimination.jl` or gray ice floe leads image in `segmentation_f.jl`
- `cloudmask`: bitmatrix cloudmask for region of interest
- `ice_labels`: vector if pixel coordinates output from `find_ice_labels.jl`

"""
function segmented_ice_cloudmasking(
    gray_image::Matrix{Gray{Float64}},
    cloudmask::BitMatrix,
    ice_labels::Union{Vector{Int64},BitMatrix,AbstractArray{<:Gray}},
)::BitMatrix
    segmented_ice = IceFloeTracker.kmeans_segmentation(gray_image, ice_labels)
    segmented_ice_cloudmasked = deepcopy(segmented_ice)
    segmented_ice_cloudmasked[cloudmask] .= 0
    return segmented_ice_cloudmasked
end

"""
    segmentation_A(segmented_ice_cloudmasked; min_opening_area)

Apply k-means segmentation to a gray image to isolate a cluster group representing sea ice. Returns an image segmented and processed as well as an intermediate files needed for downstream functions.

# Arguments

- `segmented_ice_cloudmask`: bitmatrix with open water/clouds = 0, ice = 1, output from `segmented_ice_cloudmasking()`
- `min_opening_area`: minimum size of pixels to use during morphological opening
- `fill_range`: range of values dictating the size of holes to fill

"""
function segmentation_A(
    segmented_ice_cloudmasked::BitMatrix; min_opening_area::Real=50
)::BitMatrix
    segmented_ice_opened = ImageMorphology.area_opening(
        segmented_ice_cloudmasked; min_area=min_opening_area
    )

    IceFloeTracker.hbreak!(segmented_ice_opened)

    segmented_opened_branched = IceFloeTracker.branch(segmented_ice_opened)

    segmented_bridged = IceFloeTracker.bridge(segmented_opened_branched)

    segmented_ice_filled = IceFloeTracker.fill_holes(segmented_bridged)

    diff_matrix = segmented_ice_opened .!= segmented_ice_filled

    segmented_A = segmented_ice_cloudmasked .|| diff_matrix

    return segmented_A
end

function get_holes(img, min_opening_area=20, se=IceFloeTracker.se_disk4())
    _img = ImageMorphology.area_opening(img; min_area=min_opening_area)
    IceFloeTracker.hbreak!(_img)

    out = branchbridge(_img)
    out = IceFloeTracker.opening(out, centered(se))
    out = IceFloeTracker.fill_holes(out)

    return out .!= _img
end

function fillholes!(img)
    img[get_holes(img)] .= true
    return nothing
end

function get_segment_mask(ice_mask, tiled_binmask)
    # TODO: Threads.@threads # sometimes crashes (too much memory?)
    for img in (ice_mask, tiled_binmask)
        fillholes!(img)
        img .= watershed1(img)
    end
    segment_mask = ice_mask .&& tiled_binmask
    return segment_mask
end

function branchbridge(img)
    img = IceFloeTracker.branch(img)
    img = IceFloeTracker.bridge(img)
    return img
end

"""
    imgradientmag(img)

Compute the gradient magnitude of an image using the Sobel operator.
"""
function imgradientmag(img)
    h = centered([-1 0 1; -2 0 2; -1 0 1]')
    Gx_future = Threads.@spawn IceFloeTracker.imfilter(img, h', "replicate")
    Gy_future = Threads.@spawn IceFloeTracker.imfilter(img, h, "replicate")
    Gx = fetch(Gx_future)
    Gy = fetch(Gy_future)
    return hypot.(Gx, Gy)
end