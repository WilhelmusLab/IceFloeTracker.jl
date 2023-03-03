"""
    watershed_ice_floes(intermediate_segmentation_image;)

Performs image processing and watershed segmentation with intermediate files from segmentation_b.jl or segmentation_c.jl to further isolate ice floes, returning a binary segmentation mask indicating potential sparse boundaries of ice floes.

# Arguments
-`intermediate_segmentation_image`: binary cloudmasked and landmasked intermediate file from segmentation B, either `SegB.not_ice_bit` or `SegB.ice_mask`

"""
function watershed_ice_floes(intermediate_segmentation_image::BitMatrix)::BitMatrix
    features = Images.feature_transform(.!intermediate_segmentation_image)
    distances = 1 .- Images.distance_transform(features)
    seg_mask = ImageSegmentation.hmin_transform(distances, 2)
    seg_mask_bool = seg_mask .< 1
    markers = Images.label_components(seg_mask_bool)
    segment = ImageSegmentation.watershed(distances, markers)
    labels = ImageSegmentation.labels_map(segment)
    borders = Images.isboundary(labels)
    return borders
end

"""
    segmentation_B(sharpened_image, cloudmask, segmented_a_ice_mask, struct_elem; fill_range, isolation_threshold, alpha_level, adjusted_ice_threshold)

Performs image processing and morphological filtering with intermediate files from normalization.jl and segmentation_A to further isolate ice floes, returning a mask of potential ice.

# Arguments
- `sharpened_image`: non-cloudmasked but sharpened image, output from `normalization.jl`
- `cloudmask`:  bitmatrix cloudmask for region of interest
- `segmented_a_ice_mask`: binary cloudmasked ice mask from `segmentation_a_direct.jl`
- `struct_elem`: structuring element for dilation
- `fill_range`: range of values dictating the size of holes to fill
- `isolation_threshold`: threshold used to isolated pixels from `sharpened_image`; between 0-1
- `alpha_level`: alpha threshold used to adjust contrast
- `gamma_factor`: amount of gamma adjustment 
- `adjusted_ice_threshold`: threshold used to set ice equal to one after gamma adjustment

"""
function segmentation_B(
    sharpened_image::Matrix{Gray{Float64}},
    cloudmask::BitMatrix,
    segmented_a_ice_mask::BitMatrix,
    struct_elem::ImageMorphology.MorphologySEArray{2}=strel_diamond((3, 3));
    fill_range::Tuple=(0, 1),
    isolation_threshold::Float64=0.4,
    alpha_level::Float64=0.5,
    gamma_factor::Float64=2.5,
    adjusted_ice_threshold::Float64=0.05,
)

    ## Process sharpened image
    not_ice_mask = deepcopy(sharpened_image)
    not_ice_mask[not_ice_mask .< isolation_threshold] .= 0
    not_ice_mask .= (not_ice_mask .* 0.3) .+ sharpened_image
    adjusted_sharpened = (
        (1 - alpha_level) .* sharpened_image .+ alpha_level .* not_ice_mask
    )
    gamma_adjusted_sharpened = ImageContrastAdjustment.adjust_histogram(
        adjusted_sharpened, GammaCorrection(; gamma=gamma_factor)
    )
    gamma_adjusted_sharpened_cloudmasked = IceFloeTracker.apply_cloudmask(
        gamma_adjusted_sharpened, cloudmask
    )
    segb_filled =
        .!ImageMorphology.imfill(
            gamma_adjusted_sharpened_cloudmasked .<= adjusted_ice_threshold, fill_range
        )

    ## Process ice mask
    segb_ice = MorphSE.closing(segmented_a_ice_mask, struct_elem) .* segb_filled

    ice_intersect = (segb_filled .* segb_ice)

    not_ice_bit = not_ice_mask .> 0.499 # this threshold converts pixels to zeros and ones at midpoint of the range
    segb_ice .= watershed_ice_floes(not_ice_bit)
    ice_intersect .= watershed_ice_floes(ice_intersect)

    watershed_intersect = segb_ice .* ice_intersect

    return (;
        :not_ice => map(clamp01nan, not_ice_mask)::Matrix{Gray{Float64}},
        :ice_intersect => ice_intersect::BitMatrix,
        :watershed_intersect => watershed_intersect::BitMatrix,
    )
end
