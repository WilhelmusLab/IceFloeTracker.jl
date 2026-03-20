module LopezAcosta2019

export Segment, IceDetectionLopezAcosta2019

import Images:
    Images,
    AbstractGray,
    AbstractRGB,
    adjust_histogram!,
    TransparentRGB,
    TransparentGray,
    mreconstruct!,
    mreconstruct,
    feature_transform,
    distance_transform,
    hmin_transform,
    label_components,
    watershed,
    labels_map,
    isboundary,
    SegmentedImage,
    segment_mean,
    float64,
    channelview,
    build_histogram,
    adjust_histogram,
    imfill,
    opening,
    closing,
    feature_transform,
    distance_transform,
    hmin_transform,
    clamp01nan,
    area_opening,
    area_opening!,
    dilate,
    strel_diamond,
    complement,
    bothat,
    AdaptiveEqualization,
    colorview,
    Gray,
    AbstractRGB,
    RGB,
    GammaCorrection,
    centered,
    red,
    green,
    blue

import Peaks: findmaxima
import StatsBase: kurtosis, skewness, mean, std

import ..Filtering:
    nonlinear_diffusion, PeronaMalikDiffusion, unsharp_mask, channelwise_adapthisteq
import ..Morphology: hbreak, hbreak!, branch, bridge, fill_holes, se_disk4
import ..Preprocessing:
    make_landmask_se,
    create_landmask,
    create_cloudmask,
    apply_landmask,
    apply_landmask!,
    apply_cloudmask,
    apply_cloudmask!,
    LopezAcostaCloudMask

import ..Segmentation:
    IceFloeSegmentationAlgorithm,
    find_ice_mask,
    kmeans_binarization,
    IceDetectionFirstNonZeroAlgorithm,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionThresholdMODIS721

import ..ImageUtils: imbrighten

""" 
    LopezAcosta2019.Segment(
        coastal_buffer_structuring_element::AbstractMatrix{Bool} = make_landmask_se()
        cloud_mask_algorithm = LopezAcostaCloudMask()
        diffusion_algorithm = PeronaMalikDiffusion(0.1, 0.1, 5, "exponential")
        adapthisteq_params = (
            nbins=256,
            rblocks=8, # matlab default is 8 CP
            cblocks=8, # matlab default is 8 CP
            clip=0.95,  # matlab default is 0.01 CP, which should be the same as clip=0.99
        )
        unsharp_mask_params = (smoothing_param=10, intensity=0.5)
        kmeans_params = (k=4, maxiter=50, random_seed=45)
        cluster_selection_algorithm = IceDetectionLopezAcosta2019()
        segB_params = (
            isolation_threshold=0.4,
            brightening_factor=0.3,
            gamma_factor=2.5,
            adjusted_ice_threshold=0.05,
            fill_range_max=1,
            alpha_level=0.5
        )

Segmentation algorithm for sea ice floe identification based on Lopez-Acosta 2019, 2021. The basic procedure is as follows:
1. Preprocess the image using diffusion, adaptive histogram equalization, and unsharp masking.
2. Produce a set of binary classified images using (a) global thresholds, (b) k-means cluster selection, and (c) global thresholds
   on contrast-adjusted images. Use image morphology to clean and merge individual classified images.
3. Find the shared boundaries across the segmentation methods.
4. Use the initial classification and shared boundaries to refine the preprocessed image.
5. Perform k-means binarization on the refined image
6. Use image morphology to separate ice floes.

## Arguments
- `cloud_mask_algorithm`: An `AbstractCloudMaskAlgorithm`. Defaults to [`LopezAcostaCloudMask`](@ref)
- `diffusion_algorithm`: An `AbstractDiffusionAlgorithm`. Defaults to [`PeronaMalikDiffusion`](@ref)
- `adapthisteq_params`: Parameters for the adaptive histogram AdaptiveEqualization. 
- `unsharp_mask_params`: Parameters for [`unsharp_mask`](@ref)
- `kmeans_params`: Parameters for [`kmeans_binarization`](@ref)
- `cluster_selection_algorithm`: An [`IceDetectionAlgorithm`](@ref), which takes the falsecolor image as an input and produces 
   a binary image with likely ice floe pixels set to `true`.
- `segB_params`: A collection of parameters for the second segmentation stage. `isolation_threshold` is a global threshold for
   selecting likely ice; `brightening_factor` is a percentage to increase the brightness of ice regions, `gamma_factor` is an
   input to the ImageContrastAdjustment GammaCorrection algorithm, `adjusted_ice_threshold` is a global threshold for the internal
   adjusted image, `alpha_level` controls the amount of the brightened image to use, and `fill_range_max` is the largest dark spot
   to fill in the binarized result.

## Returns
A segmented image with candidate sea ice floes.

Note: This algorithm is under active development and the API will change in a future release.
"""
@kwdef struct Segment <: IceFloeSegmentationAlgorithm
    coastal_buffer_structuring_element::AbstractMatrix{Bool} = make_landmask_se()
    cloud_mask_algorithm = LopezAcostaCloudMask()
    diffusion_algorithm = PeronaMalikDiffusion(0.1, 0.1, 5, "exponential")
    adapthisteq_params = (
        nbins=256,
        rblocks=8, # matlab default is 8 CP
        cblocks=8, # matlab default is 8 CP
        clip=0.95,  # matlab default is 0.01 CP, which should be the same as clip=0.99
    )
    unsharp_mask_params = (smoothing_param=10, intensity=0.5)
    kmeans_params = (k=4, maxiter=50, random_seed=45)
    cluster_selection_algorithm = IceDetectionLopezAcosta2019()
    segB_params = (
        isolation_threshold=0.4,
        brightening_factor=0.3,
        gamma_factor=2.5,
        adjusted_ice_threshold=0.05,
        fill_range_max=1,
        alpha_level=0.5,
    )
    segF_params = (k=3, se=strel_diamond((5, 5)), min_area_opening=20)
    floe_splitting_settings = (
        max_fill_area=1, min_area_opening=20, opening_strel=se_disk4()
    )
end

function (p::Segment)(
    truecolor::T,
    falsecolor::T,
    landmask::U;
    intermediate_results_callback::Union{Nothing,Function}=nothing,
) where {T<:AbstractMatrix{<:AbstractRGB},U<:AbstractMatrix}
    @info "building landmask and coastal buffer mask"
    _landmask_temp = create_landmask(
        float64.(landmask), p.coastal_buffer_structuring_element
    )
    landmask = _landmask_temp.non_dilated
    coastal_buffer_mask = _landmask_temp.dilated
    return p(
        truecolor,
        falsecolor,
        landmask,
        coastal_buffer_mask;
        intermediate_results_callback=intermediate_results_callback,
    )
end

function (p::Segment)(
    truecolor::T,
    falsecolor::T,
    landmask::U,
    coastal_buffer_mask::U;
    intermediate_results_callback::Union{Nothing,Function}=nothing,
) where {T<:AbstractMatrix{<:AbstractRGB},U<:BitMatrix}

    # Move these conversions down through the function as each step gets support for 
    # the full range of image formats
    truecolor_image = float64.(truecolor)
    falsecolor_image = float64.(falsecolor)

    @info "Building cloudmask"
    # TODO: Make sure tests aren't over-sensitive to roundoff errors for Float32 vs Float64
    cloudmask = create_cloudmask(falsecolor_image)

    # 2. Intermediate images
    fc_masked = apply_landmask(falsecolor_image, coastal_buffer_mask)

    @info "Preprocessing truecolor image"
    # nonlinear diffusion
    apply_landmask!(truecolor_image, landmask)
    sharpened_truecolor_image = nonlinear_diffusion(truecolor_image, p.diffusion_algorithm)

    sharpened_truecolor_image .= channelwise_adapthisteq(
        sharpened_truecolor_image;
        nbins=p.adapthisteq_params.nbins,
        rblocks=p.adapthisteq_params.rblocks,
        cblocks=p.adapthisteq_params.cblocks,
        clip=p.adapthisteq_params.clip,
    )

    sharpened_grayscale_image = unsharp_mask(
        Gray.(sharpened_truecolor_image),
        p.unsharp_mask_params.smoothing_param,
        p.unsharp_mask_params.intensity,
    )
    apply_landmask!(sharpened_grayscale_image, coastal_buffer_mask)

    # 3. Segmentation
    @info "Segmenting floes part 1/3"
    # The first segmentation routine uses "discriminate_ice_water" to enhance the contrast in the grayscale image.
    # Then, k-means clustering is used to select sea ice floes, and morphological cleanup is applied.
    ice_water_discrim = discriminate_ice_water(
        sharpened_grayscale_image, fc_masked, coastal_buffer_mask, cloudmask
    )
    segmentation_A =
        kmeans_binarization(
            ice_water_discrim,
            fc_masked;
            k=p.kmeans_params.k,
            maxiter=p.kmeans_params.maxiter,
            random_seed=p.kmeans_params.random_seed,
            cluster_selection_algorithm=p.cluster_selection_algorithm,
        ) |> clean_binary_floes

    # Potential upgrade: Remove segments of the k-means result which are all cloud. However the 
    # small isolated clouds could be filled if surrounded by a single segment.
    apply_cloudmask!(segmentation_A, cloudmask)

    @info "Segmenting floes part 2/3"
    # The second segmentation routine uses imbrighten to increase contrast between ice floes
    # and the background. It uses a simple threshold-based mask to select where to brighten.
    # Then, the segB_binarize function uses gamma correction to increase contrast before 
    # a second binary threshold is applied.

    prelim_binarized = sharpened_grayscale_image .> p.segB_params.isolation_threshold

    brightened_gray = imbrighten(
        sharpened_grayscale_image, prelim_binarized, 1 + p.segB_params.brightening_factor
    )

    segmentation_B = segB_binarize(
        sharpened_grayscale_image,
        brightened_gray,
        cloudmask;
        gamma_factor=p.segB_params.gamma_factor,
        adjusted_ice_threshold=p.segB_params.adjusted_ice_threshold,
        fill_range=(0, p.segB_params.fill_range_max),
        alpha_level=p.segB_params.alpha_level,
    )

    # Simple join of segmentations results
    ice_intersect = segmentation_A .* segmentation_B

    # Process watershed in parallel using Folds
    @info "Building watersheds"
    watersheds_segB = [
        watershed_ice_floes(prelim_binarized), watershed_ice_floes(ice_intersect)
    ]
    watersheds_product = watershed_product(watersheds_segB...)

    # segmentation_F
    # TODO: @hollandjg find out why segF is more dilated
    @info "Segmenting floes part 3/3"

    # segmentation_F
    @info "Segmenting floes part 3/3"
    morphed_grayscale = reconstruct_and_mask(
        brightened_gray,
        watersheds_product,
        ice_intersect;
        se=p.segF_params.se,
        min_area_opening=p.segF_params.min_area_opening,
    )
    # kmeans binarization, again
    segF_binarized =
        kmeans_binarization(
            morphed_grayscale,
            fc_masked;
            k=p.segF_params.k,
            cluster_selection_algorithm=p.cluster_selection_algorithm,
        ) .* .!watersheds_product

    @info "Splitting floes"
    segF = morph_split_floes(segF_binarized, cloudmask; p.floe_splitting_settings...)
    labels = label_components(segF)

    # Return the original truecolor image, segmented
    segments = SegmentedImage(truecolor, labels)

    if !isnothing(intermediate_results_callback)
        segments_truecolor = SegmentedImage(truecolor, labels)
        segments_falsecolor = SegmentedImage(falsecolor, labels)
        intermediate_results_callback(;
            truecolor,
            falsecolor,
            landmask,
            coastal_buffer_mask,
            cloudmask,
            ice_mask=IceDetectionLopezAcosta2019()(fc_masked),
            sharpened_grayscale_image=sharpened_grayscale_image,
            ice_water_discrim=ice_water_discrim,
            segA=segmentation_A,
            segB=segmentation_B,
            segAB_intersect=ice_intersect,
            watersheds_segB_product=watersheds_product,
            final_floes=segF,
            labels=labels,
            segment_mean_truecolor=map( # TODO Add "view_seg" code snippet
                i -> segment_mean(segments_truecolor, i),
                labels_map(segments_truecolor),
            ),
            segment_mean_falsecolor=map(
                i -> segment_mean(segments_falsecolor, i), labels_map(segments_falsecolor)
            ), # Add figure that overlays the segments
        )
    end
    return segments
end

"""
    discriminate_ice_water(
        sharpened_grayscale_image,
        falsecolor_image,
        landmask::T,
        cloudmask::T,
        floes_threshold::Float64=Float64(100 / 255),
        mask_clouds_lower::Float64=Float64(17 / 255),
        mask_clouds_upper::Float64=Float64(30 / 255),
        kurt_thresh_lower::Real=2,
        kurt_thresh_upper::Real=8,
        skew_thresh::Real=4,
        st_dev_thresh_lower::Float64=Float64(84 / 255),
        st_dev_thresh_upper::Float64=Float64(98.9 / 255),
        clouds_ratio_threshold::Float64=0.02,
        differ_threshold::Float64=0.6
    )

Generates an image with ice floes apparent after filtering and combining previously processed versions of falsecolor and truecolor images from the same region of interest. Returns an image ready for segmentation to isolate floes.


# Arguments
- `sharpened_grayscale_image`: Grayscale image after preprocessing
- `falsecolor_landmasked`: MODIS 7-2-1 falsecolor image after application of landmask
- `landmask`: Landmask to be used in the reconstruction function
- `cloudmask`: Cloud mask
- `floes_threshold`: Minimum band 2 and band 1 brightness for possible ice floes
- `mask_clouds_lower`: lower heuristic applied to mask out clouds
- `mask_clouds_upper`: upper heuristic applied to mask out clouds
- `kurt_thresh_lower`: lower heuristic used to set pixel value threshold based on kurtosis in histogram
- `kurt_thresh_upper`: upper heuristic used to set pixel value threshold based on kurtosis in histogram
- `skew_thresh`: heuristic used to set pixel value threshold based on skewness in histogram
- `st_dev_thresh_lower`: lower heuristic used to set pixel value threshold based on standard deviation in histogram
- `st_dev_thresh_upper`: upper heuristic used to set pixel value threshold based on standard deviation in histogram
- `clouds2_threshold`: heuristic used to set pixel value threshold based on ratio of clouds
- `differ_threshold`: heuristic used to calculate proportional intensity in histogram
"""
function discriminate_ice_water(
    sharpened_grayscale_image, #::AbstractArray{AbstractGray}, #dmw: discrim-ice-water test fails here
    landmasked_falsecolor_image, #::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    landmask::T,
    cloudmask::T,
    floes_threshold::Float64=Float64(100 / 255),
    mask_clouds_lower::Float64=Float64(17 / 255),
    mask_clouds_upper::Float64=Float64(30 / 255),
    kurt_thresh_lower::Real=2,
    kurt_thresh_upper::Real=8,
    skew_thresh::Real=4,
    st_dev_thresh_lower::Float64=Float64(84 / 255),
    st_dev_thresh_upper::Float64=Float64(98.9 / 255),
    clouds_ratio_threshold::Float64=0.02,
    differ_threshold::Float64=0.6,
)::AbstractMatrix where {T<:AbstractArray{Bool}}

    # First step: Grayscale reconstruction, creating an inverted and smoothed image.
    fc_landmasked = landmasked_falsecolor_image # shorten name for convenience
    morphed_grayscale = _reconstruct(sharpened_grayscale_image, landmask)

    # This next section is all here to find a threshold value for masking. Only three levels are selected,
    # which makes me think that we'd do better to use a percentile function or otherwise continuous estimate
    # of a threshold from the data, rather than 3 fixed steps.
    b7_landmasked = Gray.(red.(fc_landmasked)) # formerly falsecolor_image_band7 -> image_cloudless (but it does have clouds???)
    b7_landmasked_cloudmasked = apply_cloudmask(b7_landmasked, cloudmask) # formerly clouds_channel -> image_clouds

    # Select pixels greater than intensity 100 in bands 2 and 1
    b2_landmasked = green.(fc_landmasked)
    b1_landmasked = blue.(fc_landmasked)
    b2_subset = b2_landmasked[b2_landmasked .> floes_threshold]
    b1_subset = b1_landmasked[b1_landmasked .> floes_threshold]

    # if b2_subset, b1_subset are empty, then there are no floes to find and we could end the algorithm right there.
    # question is whether we should return a blank image, or if we should return the unmodified grayscale image.
    # Also, in this case, is blank an image of 1s or an image of 0s?
    length(b2_subset) < 10 && return morphed_grayscale

    # Compute "proportional intensity", a measure of the prominence of a peak
    # The nbins value is just the number of gray levels larger than the floes threshold.
    nbins = round(Int64, 255 * (1 - floes_threshold))
    _, floes_bin_counts = build_histogram(b2_subset, nbins)
    _, vals = findmaxima(floes_bin_counts)
    differ = vals / (maximum(vals))
    proportional_intensity = sum(differ .> differ_threshold) / length(differ)

    # compute kurtosis, skewness, and standard deviation to use in threshold filtering
    kurt_band_2 = kurtosis(b2_subset)
    skew_band_2 = skewness(b2_subset)
    kurt_band_1 = kurtosis(b1_subset)
    standard_dev = std(vec(morphed_grayscale))

    # The clouds ratio was computed on the whole area, which means that 
    # there will be errors near the land mask. Correcting this may make it 
    # have different results than the Matlab version.
    clouds_ratio = mean(b7_landmasked_cloudmasked[.!landmask] .> 0)

    # It may be worthwhile to take a random sample of scenes and test what the kurtosis, skew, and intensity are.
    # These values are likely to vary with the size of the image. Both band 1 and band 2 are used, though they 
    # are highly correlated with each other.
    threshold_50_check = _check_threshold_50(
        kurt_band_1,
        kurt_band_2,
        kurt_thresh_lower,
        kurt_thresh_upper,
        skew_band_2,
        skew_thresh,
        proportional_intensity,
    )

    # This method uses standard deviation of the grayscale image instead of the band 1 / band 2 values.
    threshold_130_check = _check_threshold_130(
        clouds_ratio,
        clouds_ratio_threshold,
        standard_dev,
        st_dev_thresh_lower,
        st_dev_thresh_upper,
    )

    # If neither check passes, then a middle value is used (not exactly the middle, but close).
    if threshold_50_check
        THRESH = 50 / 255
    elseif threshold_130_check
        THRESH = 130 / 255
    else
        THRESH = 80 / 255 #intensity value of 80
    end

    # Values less than the threshold are set to 0. Essentially, it's flattening out the "ice" portion of the image.
    morphed_image_copy = copy(morphed_grayscale)
    morphed_image_copy[morphed_grayscale .<= THRESH] .= 0

    # Finally, there is a mask applied to a narrow band of band7 values. Removing
    # the next section *does not affect the tests* which makes me wonder if we need it.
    # What's happening here is that the cloud_threshold variable is selecting anything darker than
    # clouds lower, or brighter than clouds lower.
    # So the b7_landmasked getting multiplied by not(cloud_threshold) means that we're selecting all
    # the intermediate values: anything between mask_clouds_lower and mask_clouds_upper.
    # This is essentially speckle noise, and has minimal effect on the results.

    _cloud_threshold = (
        b7_landmasked_cloudmasked .< mask_clouds_lower .||
        b7_landmasked_cloudmasked .> mask_clouds_upper
    )

    # reusing image_cloudless - used to be band7_masked
    # I think the names were backwards: image_cloudless had not been cloudmasked, and image_clouds had been.
    # So the question is if they're swapped in 305 and 311 also.
    @. b7_landmasked = b7_landmasked * !_cloud_threshold

    # Check to see if selecting indices to set to 0 would be equivalent
    @. morphed_image_copy = clamp01nan(morphed_image_copy - (b7_landmasked * 3))

    return morphed_image_copy
end

function _check_threshold_50(
    kurt_band_1,
    kurt_band_2,
    kurt_thresh_lower,
    kurt_thresh_upper,
    skew_band_2,
    skew_thresh,
    proportional_intensity,
)
    return ( # intensity value of 50
        (
            (kurt_band_2 > kurt_thresh_upper) ||
            (kurt_band_2 < kurt_thresh_lower) && (kurt_band_1 > kurt_thresh_upper)
        ) ||
        (
            (kurt_band_2 < kurt_thresh_lower) &&
            (skew_band_2 < skew_thresh) &&
            proportional_intensity < 0.1
        ) ||
        proportional_intensity < 0.01
    )
end

function _check_threshold_130(
    clouds_ratio,
    clouds_ratio_threshold,
    standard_dev,
    st_dev_thresh_lower,
    st_dev_thresh_upper,
)
    return (clouds_ratio .< clouds_ratio_threshold && standard_dev > st_dev_thresh_lower) ||
           (standard_dev > st_dev_thresh_upper)
end

"""_reconstruct(sharpened_grayscale_image, dilated_mask; strel)

Convenience function for reconstruction by dilation using the complement
of an image. Markers are computed by dilating the input image by the 
structuring element `strel` and taking the complement. The dilated landmask
is applied at the end to prevent bright regions from bleeding into the land mask.
Defaults to using a radius 5 diamond mask.
"""
function _reconstruct(sharpened_grayscale_image, dilated_mask; strel=strel_diamond((5, 5)))
    markers = complement.(dilate(sharpened_grayscale_image, strel))
    mask = complement.(sharpened_grayscale_image)
    reconstructed_grayscale = mreconstruct(dilate, markers, mask)
    apply_landmask!(reconstructed_grayscale, dilated_mask)
    return reconstructed_grayscale
end

"""
    clean_binary_floes(bw_img; min_opening_area=50)

Refine a binarized ice floe image (floes=white, leads/background/water=black) using morphological operations.

"""
function clean_binary_floes(bw_img; min_opening_area=50, se=strel_diamond((5, 5)))
    img_opened = area_opening(bw_img; min_area=min_opening_area) |> hbreak
    img_filled = branch(img_opened) |> bridge |> fill_holes
    diff_matrix = img_opened .!= img_filled
    return closing(bw_img .|| diff_matrix, se) # Original has closing here, but opening works better!
end

"""
 segB_binarize(sharpened_image, brightened_image, cloudmask;
     gamma_factor=2.5, adjusted_ice_threshold=0.05, fill_range=(0, 1), alpha_level=0.5)

Binarize the sharpened image by selective brightening, gamma correction, threshold application, and 
clean up with image hole filling.

"""
function segB_binarize(
    sharpened_image,
    brightened_image,
    cloudmask;
    gamma_factor=2.5,
    adjusted_ice_threshold=0.05,
    fill_range=(0, 1),
    alpha_level=0.5,
)
    # Weighted average between brightened image and sharpened grayscale
    adjusted_sharpened =
        (1 - alpha_level) .* sharpened_image .+ alpha_level .* brightened_image

    # Gamma correction and cloud masking
    adjust_histogram!(adjusted_sharpened, GammaCorrection(; gamma=gamma_factor))
    apply_cloudmask!(adjusted_sharpened, cloudmask) # Should this be happening here?

    # Thresholding and filling small holes (based on fill_range=(min, max))
    segB = adjusted_sharpened .<= adjusted_ice_threshold
    segb_filled = .!imfill(segB, fill_range)
    return segb_filled
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
    gray_image, falsecolor_image, cloudmask::BitMatrix
)::BitMatrix
    segmented_ice = kmeans_binarization(
        gray_image,
        falsecolor_image;
        cluster_selection_algorithm=IceDetectionLopezAcosta2019(),
    )
    segmented_ice_cloudmasked = deepcopy(segmented_ice)
    segmented_ice_cloudmasked[cloudmask] .= 0
    return segmented_ice_cloudmasked
end
"""
    watershed_ice_floes(intermediate_segmentation_image;)
Performs image processing and watershed segmentation with intermediate files from segmentation_b.jl to further isolate ice floes, returning a binary segmentation mask indicating potential sparse boundaries of ice floes.
# Arguments
-`intermediate_segmentation_image`: binary cloudmasked and landmasked intermediate file from segmentation B, either `SegB.not_ice_bit` or `SegB.ice_intersect`
"""
function watershed_ice_floes(intermediate_segmentation_image::BitMatrix)::BitMatrix
    features = feature_transform(.!intermediate_segmentation_image)
    distances = 1 .- distance_transform(features)
    seg_mask = hmin_transform(distances, 2)
    seg_mask_bool = seg_mask .> 0
    markers = label_components(seg_mask_bool)
    segment = watershed(distances, markers)
    labels = labels_map(segment)
    borders = isboundary(labels)
    return borders
end

"""
    watershed_product(watershed_B_ice_intersect, watershed_B_not_ice;)
Intersects the outputs of watershed segmentation on intermediate files from segmentation B, indicating potential sparse boundaries of ice floes.
# Arguments
- `watershed_B_ice_intersect`: binary segmentation mask from `watershed_ice_floes`
- `watershed_B_not_ice`: binary segmentation mask from `watershed_ice_floes`
"""
function watershed_product(
    watershed_B_ice_intersect::BitMatrix, watershed_B_not_ice::BitMatrix;
)::BitMatrix

    ## Intersect the two watershed files
    watershed_intersect = watershed_B_ice_intersect .* watershed_B_not_ice
    return watershed_intersect
end

"""
    _adjust_histogram(masked_view, nbins, rblocks, cblocks, clip)

Perform adaptive histogram equalization to a masked image. Wrapper for the
JuliaImages `adjust_histogram` function with `AdaptiveEqualization`, setting the minval
and maxval to the image maximum and minimum.

"""
function _adjust_histogram(masked_view, nbins, rblocks, cblocks, clip)
    return adjust_histogram(
        masked_view,
        AdaptiveEqualization(;
            nbins=nbins,
            rblocks=rblocks,
            cblocks=cblocks,
            minval=minimum(masked_view), # Could this be causing the unnatural coloration in dark image regions?
            maxval=maximum(masked_view),
            clip=clip,
        ),
    )
end

"""IceDetectionLopezAcosta2019

Application of the IceDetectionFirstNonZeroAlgorithm using two passes of 
the IceDetectionThresholdMODIS721 and one application of the IceDetectionBrightnessPeaksMODIS721.
""" # TODO: This works in the kmeans binarization but not by itself in the example notebook.
function IceDetectionLopezAcosta2019(;
    band_7_max::Float64=Float64(5 / 255),
    band_2_min::Float64=Float64(230 / 255),
    band_1_min::Float64=Float64(240 / 255),
    band_7_max_relaxed::Float64=Float64(10 / 255),
    band_1_min_relaxed::Float64=Float64(190 / 255),
    possible_ice_threshold::Float64=Float64(75 / 255),
)
    return IceDetectionFirstNonZeroAlgorithm(
        [
            IceDetectionThresholdMODIS721(;
                band_7_max=band_7_max, band_2_min=band_2_min, band_1_min=band_1_min
            ),
            IceDetectionThresholdMODIS721(;
                band_7_max=band_7_max_relaxed,
                band_2_min=band_2_min,
                band_1_min=band_1_min_relaxed,
            ),
            IceDetectionBrightnessPeaksMODIS721(;
                band_7_max=band_7_max, possible_ice_threshold=possible_ice_threshold
            ),
        ],
        10,
    )
end

"""
    reconstruct_and_mask(grayscale_img, watershed_boundary, ice_intersect;
    se=strel_diamond((5,5)), min_area_opening=20
    )

Enhance the visibility of distinct floes in the grayscale image by using grayscale reconstruction.
Updates the sea ice mask by intersecting the `ice_intersect` and the `watershed_boundary` and using 
morphological area opening.
"""
function reconstruct_and_mask(
    grayscale_img::AbstractArray{<:Union{AbstractGray,TransparentGray}},
    watershed_boundary::BitMatrix,
    ice_intersect::BitMatrix;
    se=strel_diamond((5, 5)),
    min_area_opening=20,
)
    ice_mask = .!watershed_boundary .* ice_intersect
    ice_mask .= .!area_opening(ice_mask; min_area=min_area_opening, connectivity=2)

    reconst_gray = dilate(grayscale_img, se)
    mreconstruct!(
        dilate, reconst_gray, complement.(reconst_gray), complement.(grayscale_img)
    )

    apply_landmask!(reconst_gray, ice_mask .== 0)
    return reconst_gray
end

"""
    morph_split_floes(binary_img; max_fill_area=1, min_area_opening=20, opening_strel=se_disk(4))

Separate floes in a binary image using hbreak, branch, imfill, bottom hat transform, and area opening.
Based on Lopez-Acosta et al. 2019, 2021.

"""
function morph_split_floes(
    binary_img, cloudmask; max_fill_area=1, min_area_opening=20, opening_strel=se_disk4()
)
    leads_branched = hbreak(binary_img) |> branch
    leads_filled = .!imfill(.!leads_branched, 0:max_fill_area)
    leads_opened = branch(
        area_opening(leads_filled; min_area=min_area_opening, connectivity=2)
    )

    leads_bothat = bothat(leads_opened, strel_diamond((5, 5))) .> 0
    leads = convert(BitMatrix, (complement.(leads_bothat) .* leads_opened))
    area_opening!(leads, leads; min_area=min_area_opening, connectivity=2)

    floes = (fill_holes(leads) .* .!cloudmask) |> branch
    floes_opened = opening(floes, opening_strel)
    mreconstruct!(dilate, floes_opened, floes, floes_opened)
    return floes_opened
end

end
