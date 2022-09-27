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
    struct_elem::Matrix{Bool};
    fill_range::Tuple=(0, 40),
    isolation_threshold::Float64=0.4,
    alpha_level::Float64=0.5,
    gamma_factor::Float64=2.5,
    adjusted_ice_threshold::Float64=0.2,
)::BitMatrix
    ## Process sharpened image
    not_ice_mask = .!(sharpened_image .< isolation_threshold)
    adjusted_sharpened = (1 - alpha_level) .* sharpened_image .+ alpha_level .* not_ice_mask

    gamma_adjusted_sharpened = ImageContrastAdjustment.adjust_histogram(
        adjusted_sharpened, GammaCorrection(; gamma=gamma_factor)
    )
    gamma_adjusted_sharpened_cloudmasked, _ = IceFloeTracker.apply_cloudmask(
        gamma_adjusted_sharpened, cloudmask
    )
    gamma_adjusted_sharpened_cloudmasked[gamma_adjusted_sharpened_cloudmasked .<= adjusted_ice_threshold] .=
        0
    gamma_adjusted_sharpened_cloudmasked[gamma_adjusted_sharpened_cloudmasked .> adjusted_ice_threshold] .=
        1
    gamma_adjusted_sharpened_cloudmasked_bit = convert(
        BitMatrix, gamma_adjusted_sharpened_cloudmasked
    )

    adjusted_bitmatrix = .!(gamma_adjusted_sharpened_cloudmasked_bit)

    segb_filled = ImageMorphology.imfill(adjusted_bitmatrix, fill_range)
    segb_filled = .!(segb_filled)

    ## Process ice mask
    segmented_a_ice_mask_holes = ImageMorphology.imfill(.!segmented_a_ice_mask, fill_range)
    segmented_a_ice_masked_filled = .!segmented_a_ice_mask_holes
    segb_closed = ImageMorphology.closing(segmented_a_ice_masked_filled, struct_elem)

    ## Create mask from intersect of processed images
    segb_filled_ice = (segb_filled .> 0)
    segb_closed_ice = (segb_closed .> 0)
    segmented_b_ice_intersect = (segb_filled_ice .* segb_closed_ice)

    return segmented_b_ice_intersect
end
