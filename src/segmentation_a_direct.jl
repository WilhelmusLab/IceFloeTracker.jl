"""
    remove_landmask(landmask, ice_mask)

Find the pixel indexes that are floating ice rather than soft or land ice. Returns an array of pixel indexes. 

# Arguments
- `landmask`: bitmatrix landmask for region of interest
- `ice_mask`: bitmatrix with ones equal to ice, zeros otherwise

"""
function remove_landmask(landmask::BitMatrix, ice_mask::BitMatrix)::Array{Int64}
    indexes_no_landmask = []
    land = IceFloeTracker.apply_landmask(ice_mask, landmask)
    for (idx, val) in enumerate(land)
        if val != 0
            push!(indexes_no_landmask, idx)
        end
    end
    return indexes_no_landmask
end

"""
    segmentation_A(gray_image, cloudmask, ice_labels; min_opening_area, fill_range)

Apply k-means segmentation to a gray image to isolate a cluster group representing sea ice. Returns an image segmented and processed as well as intermediate files needed for downstream functions.

# Arguments

- `gray_image`: output image from `ice-water-discrimination.jl` or gray ice floe leads image in `segmentation_f.jl`
- `cloudmask`: bitmatrix cloudmask for region of interest
- `ice_labels`: vector if pixel coordinates output from `find_ice_labels.jl`
- `min_opening_area`: minimum size of pixels to use during morphoilogical opening
- `fill_range`: range of values dictating the size of holes to fill

"""
function segmentation_A(
    gray_image::Matrix{Gray{Float64}},
    cloudmask::BitMatrix,
    ice_labels::Vector{Int64};
    min_opening_area::Real=50,
    fill_range::Tuple=(0, 50),
)::Tuple{BitMatrix,BitMatrix,BitMatrix}
    gray_image = float64.(gray_image)
    gray_image_height, gray_image_width = size(gray_image)
    gray_image_1d = reshape(gray_image, 1, gray_image_height * gray_image_width)
    println("Done with reshape")
    feature_classes = Clustering.kmeans(
        gray_image_1d, 4; maxiter=50, display=:iter, init=:kmpp
    )
    class_assignments = assignments(feature_classes)

    ## NOTE(tjd): this reshapes column major vector of kmeans classes back into original image shape
    segmented = reshape(class_assignments, gray_image_height, gray_image_width)

    nlabel = StatsBase.mode(segmented[ice_labels])

    ## Isolate ice floes and contrast from background
    segmented_ice = segmented .== nlabel
    segmented_ice_cloudmasked = segmented_ice .* cloudmask

    segmented_ice_opened = ImageMorphology.area_opening(
        segmented_ice_cloudmasked; min_area=min_opening_area
    ) #BW_test in matlab code

    segmented_ice_opened_hbreak = IceFloeTracker.hbreak!(segmented_ice_opened)

    segmented_opened_flipped = .!segmented_ice_opened_hbreak

    segmented_ice_filled = ImageMorphology.imfill(
        convert(BitMatrix, (segmented_opened_flipped)), fill_range
    ) #BW_test3 in matlab code
    println("Done filling segmented_ice")
    segmented_ice_filled_comp = complement.(segmented_ice_filled)

    diff_matrix = segmented_ice_opened .!= segmented_ice_filled_comp

    segmented_A = segmented_ice_cloudmasked .|| diff_matrix

    return segmented_ice, segmented_ice_cloudmasked, segmented_A
end
