"""
    kmeans_segmentation(gray_image, ice_labels;)

Apply k-means segmentation to a gray image to isolate a cluster group representing sea ice. Returns a binary image with ice segmented from background.

# Arguments
- `gray_image`: output image from `ice-water-discrimination.jl` or gray ice floe leads image in `segmentation_f.jl`
- `ice_labels`: vector if pixel coordinates output from `find_ice_labels.jl`

"""
function kmeans_segmentation(
    gray_image::Matrix{Gray{Float64}},
    ice_labels::Vector{Int64},
    k::Int64=4,
    maxiter::Int64=50,
)::BitMatrix
    segmented = kmeans_segmentation(gray_image; k=k, maxiter=maxiter)

    ## Isolate ice floes and contrast from background
    segmented_ice = get_segmented_ice(segmented, ice_labels)
    return segmented_ice
end

function kmeans_segmentation(
    gray_image::Matrix{Gray{Float64}}; k::Int64=4, maxiter::Int64=50
)
    Random.seed!(45)

    ## NOTE(tjd): this clusters into k classes and solves iteratively with a max of maxiter iterations
    feature_classes = Clustering.kmeans(
        vec(gray_image), k; maxiter=maxiter, display=:none, init=:kmpp
    )

    class_assignments = assignments(feature_classes)

    ## NOTE(tjd): this clusters into 4 classes and solves iteratively with a max of 50 iterations
    segmented = reshape(class_assignments, size(gray_image))

    return segmented
end

function get_segmented_ice(segmented::Matrix{Int64}, ice_labels::Vector{Int64})
    ## Isolate ice floes and contrast from background
    nlabel = StatsBase.mode(segmented[ice_labels])
    segmented_ice = segmented .== nlabel
    return segmented_ice
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
    gray_image::Matrix{Gray{Float64}}, cloudmask::BitMatrix, ice_labels::Vector{Int64}
)::BitMatrix
    segmented_ice = IceFloeTracker.kmeans_segmentation(gray_image, ice_labels)
    segmented_ice_cloudmasked = segmented_ice .* cloudmask
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

    segmented_ice_filled = IceFloeTracker.MorphSE.fill_holes(segmented_bridged)

    diff_matrix = segmented_ice_opened .!= segmented_ice_filled

    segmented_A = segmented_ice_cloudmasked .|| diff_matrix

    return segmented_A
end

function get_holes(img, min_opening_area=20, se=IceFloeTracker.se_disk4())
    img .= ImageMorphology.area_opening(img; min_area=min_opening_area)
    IceFloeTracker.hbreak!(img)

    out = branchbridge(img)
    out = IceFloeTracker.MorphSE.opening(out, centered(se))
    out = IceFloeTracker.MorphSE.fill_holes(out)

    return out .!= img
end

function fillholes!(img)
    img[get_holes(img)] .= true
    return nothing
end

function get_segment_mask(ice_mask, tiled_binmask)
    Threads.@threads for img in (ice_mask, tiled_binmask)
        fillholes!(img)
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
    Gx = imfilter(img, h', "replicate")
    Gy = imfilter(img, h, "replicate")
    return hypot.(Gx, Gy)
end
