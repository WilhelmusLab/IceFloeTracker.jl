"""
    kmeans_segmentation(gray_image, ice_labels;)

Apply k-means segmentation to a gray image to isolate a cluster group representing sea ice. Returns a binary image with ice segmented from background.

# Arguments
- `gray_image`: output image from `ice-water-discrimination.jl` or gray ice floe leads image in `segmentation_f.jl`
- `ice_labels`: vector if pixel coordinates output from `find_ice_labels.jl`

"""
function kmeans_segmentation(
    gray_image::Matrix{Gray{Float64}}, ice_labels::Vector{Int64}
)::BitMatrix
    Random.seed!(45) # this seed generates consistent clusters for the final output
    gray_image_height, gray_image_width = size(gray_image)
    gray_image_1d = vec(gray_image)
    @info("Done with reshape")

    ## NOTE(tjd): this clusters into 4 classes and solves iteratively with a max of 50 iterations
    feature_classes = Clustering.kmeans(
        gray_image_1d, 4; maxiter=50, display=:none, init=:kmpp
    )
    class_assignments = assignments(feature_classes)

    ## NOTE(tjd): this reshapes column major vector of kmeans classes back into original image shape
    segmented = reshape(class_assignments, gray_image_height, gray_image_width)

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

    #segmented_ice_filled = .!bwareamaxfilt(.!segmented_bridged)
    segmented_ice_filled = IceFloeTracker.MorphSE.fill_holes(segmented_bridged)
    @info "Done filling segmented_ice"

    diff_matrix = segmented_ice_opened .!= segmented_ice_filled

    segmented_A = segmented_ice_cloudmasked .|| diff_matrix

    return segmented_A
end

function get_holes(img, min_opening_area, se=IceFloeTracker.MorphSE.se_disk4())
    img .= ImageMorphology.area_opening(img; min_area=min_opening_area)
    IceFloeTracker.hbreak!(img)

    out = branchbridge(img)
    out = IceFloeTracker.MorphSE.opening(out, centered(se))
    out = IceFloeTracker.MorphSE.fill_holes(out)

    return out .!= img

end

function branchbridge(img)
    img = IceFloeTracker.branch(img)
    img = IceFloeTracker.bridge(img)
    return img
end

