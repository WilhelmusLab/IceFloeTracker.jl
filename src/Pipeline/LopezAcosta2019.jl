"""
Sea ice floe segmentation algorithm version 1
The MATLAB version of this algorithm was developed as an extension of 
Lopez-Acosta et al. 2019 for Dr. Rosalinda Lopez-Acosta's doctoral research.
It was used in the production of the IFT Fram Strait Dataset (Lopez-Acosta et al. 2024).
The workflow here reproduces that result to the extent possible, and allows adaptation
for differing parameter choices.
"""
module LopezAcosta2019

import Images:
    Images,
    AbstractGray,
    AbstractRGB,
    TransparentRGB,
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
    adjust_histogram!,
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
    erode,
    strel_diamond,
    complement,
    bothat,
    AdaptiveEqualization,
    colorview,
    Gray,
    AbstractRGB,
    RGB,
    GammaCorrection,
    centered

import Peaks: findmaxima
import StatsBase: kurtosis, skewness
import ..Filtering: nonlinear_diffusion, PeronaMalikDiffusion, unsharp_mask, _channelwise_adapthisteq
import ..Morphology: hbreak, hbreak!, branch, bridge, fill_holes, se_disk4
import ..Preprocessing:
    make_landmask_se,
    create_landmask,
    create_cloudmask,
    LopezAcostaCloudMask,
    create_clouds_channel,
    apply_landmask,
    apply_landmask!,
    apply_cloudmask,
    apply_cloudmask!
import ..Segmentation:
    IceFloeSegmentationAlgorithm, 
    IceDetectionThresholdMODIS721,
    IceDetectionFirstNonZeroAlgorithm,
    IceDetectionBrightnessPeaksMODIS721,
    find_ice_mask, 
    kmeans_segmentation, 
    kmeans_binarization

"""
Sample input parameters expected by the main function
"""
cloud_mask_thresholds = (
    prelim_threshold=110.0 / 255.0,
    band_7_threshold=200.0 / 255.0,
    band_2_threshold=190.0 / 255.0,
    ratio_lower=0.0,
    ratio_offset=0.0,
    ratio_upper=0.75,
)

diffusion_parameters = (lambda=0.1, kappa=0.1, niters=5, g="exponential")

ice_masks_params = (
    band_7_max=5/255,
    band_2_min=230/255,
    band_1_min=240/255,
    band_7_max_relaxed=10 / 255,
    band_1_min_relaxed=190 / 255,
    possible_ice_threshold=75 / 255,
)

# Ice label algorithm -- identifies bright ice pixels for use in selecting the 
# k-means cluster which corresponds to the brightest floes.

function IceDetectionLopezAcosta2019(;
    band_7_max::Float64=Float64(5 / 255),
    band_2_min::Float64=Float64(230 / 255),
    band_1_min::Float64=Float64(240 / 255),
    band_7_max_relaxed::Float64=Float64(10 / 255),
    band_1_min_relaxed::Float64=Float64(190 / 255),
    possible_ice_threshold::Float64=Float64(75 / 255),
)
    return IceDetectionFirstNonZeroAlgorithm([
        IceDetectionThresholdMODIS721(;
            band_7_max=band_7_max,
            band_2_min=band_2_min,
            band_1_min=band_1_min
        ),
        IceDetectionThresholdMODIS721(;
            band_7_max=band_7_max_relaxed,
            band_2_min=band_2_min,
            band_1_min=band_1_min_relaxed,
        ),
        IceDetectionBrightnessPeaksMODIS721(;
            band_7_max=band_7_max,
            possible_ice_threshold=possible_ice_threshold
        ),
    ])
end

@kwdef struct Segment <: IceFloeSegmentationAlgorithm
    coastal_buffer_structuring_element::AbstractMatrix{Bool} = make_landmask_se()
    cloud_mask_algorithm = LopezAcostaCloudMask(cloud_mask_thresholds...)
    diffusion_algorithm = PeronaMalikDiffusion(diffusion_parameters...)
    adapthisteq_params = (
        nbins=256,
        rblocks=8, # matlab default is 8 CP
        cblocks=8, # matlab default is 8 CP
        clip=0.9,  # matlab default is 0.01 CP, which should be the same as clip=0.99
    )
    unsharp_mask_params = (smoothing_param=10, intensity=2)
    reconstruct_strel = strel_diamond((5, 5)) # Structuring element used for enhancing foreground/background contrast
    kmeans_params = (k=4, maxiter=50, random_seed=45)
    cluster_selection_algorithm = IceDetectionLopezAcosta2019()
    segmentation_b_params = (isolation_threshold=0.4, struct_elem=strel_diamond((3, 3)))
end

# function Segment(; landmask_structuring_element=make_landmask_se())
#     return Segment(landmask_structuring_element)
# end

# dmw: how can we include a pre-generated landmask? if we use map, does the compiler recognize that it's a repeat computation?
# one option would be to have args landmask_dilated and landmask_nondilated, and a method that would create_landmask if only one landmask
# is inputted. 
function (p::Segment)(
    truecolor::T,
    falsecolor::T,
    landmask::U;
    intermediate_results_callback::Union{Nothing,Function}=nothing,
) where {T<:AbstractMatrix{<:AbstractRGB},U<:AbstractMatrix}
    @info "building landmask and coastal buffer mask"
    landmask, coastal_buffer_mask = create_landmask(
        float64.(landmask), p.coastal_buffer_structuring_element
    )
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
    
    landmask_image = float64.(landmask)

    @info "building landmask"
    landmask_imgs = create_landmask(landmask_image, p.landmask_structuring_element)

    @info "Building cloudmask"
    # TODO: @hollandjg track down why the cloudmask is different for float32 vs float64 input images
    cloudmask = create_cloudmask(falsecolor_image, p.cloud_mask_algorithm)

    fc_landmasked = apply_landmask(falsecolor_image, landmask_imgs.dilated)

    @info "Preprocessing truecolor image"
    # nonlinear diffusion
    apply_landmask!(truecolor_image, landmask_imgs.non_dilated)
    sharpened_truecolor_image = nonlinear_diffusion(
        truecolor_image, p.diffusion_algorithm
    )

    # q: do we need to keep the sharpened truecolor image? or just the grayscale?
    sharpened_truecolor_image .= _channelwise_adapthisteq(sharpened_truecolor_image;
        nbins=p.adapthisteq_params.nbins,
        rblocks=p.adapthisteq_params.rblocks,
        cblocks=p.adapthisteq_params.cblocks,
        clip=p.adapthisteq_params.clip
    )

    sharpened_grayscale_image = unsharp_mask(
        Gray.(sharpened_truecolor_image),
        p.unsharp_mask_params.smoothing_param,
        p.unsharp_mask_params.intensity,
    )
    apply_landmask!(sharpened_grayscale_image, landmask_imgs.dilated)

    @info "Segmentation method 1: Grayscale reconstruction + k-means binarization"
    # reconstruction of the complement, producing image with dark floes and bright leads
    # dmw: do any of the in-between steps get reused? 
    # separating them like this lets us export individual steps, but it's more verbose.
    markers = complement.(dilate(sharpened_grayscale_image, p.reconstruct_strel))
    mask = complement.(sharpened_grayscale_image)
    reconstructed_grayscale = mreconstruct(dilate, markers, mask)
    apply_landmask!(reconstructed_grayscale, landmask_imgs.dilated)

    # further enhancing darkness of floes and brightness of leads.
    ice_water_discrim = discriminate_ice_water(
        falsecolor_image, reconstructed_grayscale, landmask_imgs.dilated, cloudmask
    )

    kmeans_result = kmeans_binarization(
        ice_water_discrim,
        fc_landmasked;
        k=p.kmeans_params.k,
        maxiter=p.kmeans_params.maxiter,
        random_seed=p.kmeans_params.random_seed,
        cluster_selection_algorithm=p.cluster_selection_algorithm
        ) |> clean_binary_floes

    # check: are there any regions that are nonzero under the cloudmask, since it was applied in discriminate ice water?
    segA = apply_cloudmask(kmeans_result, cloudmask) 

    @info "Segmentation method 2: Thresholding brightened image"

    # Brighten ice and set areas darker than a threshold to 0
    # General purpose with brighten.jl?
    # TODO: determine name for parameter, operation for brightening
    threshold_mask = sharpened_grayscale_image .> p.segmentation_b_params.isolation_threshold
    brightened_image = (sharpened_grayscale_image .* 1.3) .* threshold_mask
    clamp!(brightened_image, 0, 1)
  
    @info "Segmentation method 3: Thresholding brightened image after gamma correction"
    # The brightening could happen inside the algorithm, since it is so simple.
    segB = segB_binarize(sharpened_grayscale_image, brightened_image, cloudmask)
    
    @info "Merging segmentation results"
    segAB_intersect = closing(segA, p.segmentation_b_params.struct_elem) .* segB

    # Julia's watershed boundaries are larger than in matlab, which may be an issue. Can 
    # we improve on the boundary identification by only drawing the boundary on the non-background regions?
    watersheds_product = watershed_ice_floes(threshold_mask) .* watershed_ice_floes(segAB_intersect)

    # segmentation_F
    @info "Segmenting floes part 3/3"
    # rename
    segF = segmentation_F(
        brightened_image,
        segAB_intersect, 
        watersheds_product,
        fc_landmasked,
        cloudmask, # potentially unused
        landmask_imgs.dilated; # potentially unused
        cluster_selection_algorithm=p.cluster_selection_algorithm
    )

    binary_cleaned = morphological_cleanup_final(
        segF,
        cloudmask;
        fill_range=(0,1),
        min_area_opening=20,
        bothat_strel=strel_diamond((5,5)),
        opening_strel=centered(se_disk4()), # consider replacement with optimized strel 
    )

    @info "Labeling floes"
    labels = label_components(binary_cleaned)

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
            cloudmask=cloudmask,
            ice_mask=IceDetectionLopezAcosta2019()(fc_masked),
            sharpened_truecolor_image=sharpened_truecolor_image,
            sharpened_gray_truecolor_image=sharpened_grayscale_image,
            normalized_image=reconstructed_grayscale,
            ice_water_discrim=ice_water_discrim,
            segA=segA,
            segB=segB,
            watersheds_segB_product=watersheds_product,
            segF=segF,
            labels=labels,
            segment_mean_truecolor=map(
                i -> segment_mean(segments_truecolor, i), labels_map(segments_truecolor)
            ),
            segment_mean_falsecolor=map(
                i -> segment_mean(segments_falsecolor, i), labels_map(segments_falsecolor)
            ),
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
- `nbins`: number of bins during histogram build

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
    differ_threshold::Float64=0.6
)::AbstractMatrix where {T<:AbstractArray{Bool}} #dmw extend to allow BitMatrix or Bool if that throws an error

    # First step: Grayscale reconstruction, creating an inverted and smoothed image.
    fc_landmasked = landmasked_falsecolor_image # shorten name for convenience
    morphed_grayscale = _reconstruct(sharpened_grayscale_image, landmask)

    # Second step: Find a threshold value to mask. There are only three levels considered,
    # which makes me think that we'd do better to use a percentile function or otherwise 
    # continuous estimate of a threshold from the data, rather than 3 fixed steps.
    # The language is confusing because image_clouds / image_cloudless are backwards: image_clouds has had the image
    # masked by the cloud mask, while image_cloudless has not been cloudmasked.

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
    nbins = round(Int64, 255*(1 - floes_threshold))
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
    # clouds lower, or brighter than clouds lower. So the b7_landmasked getting multiplied by 
    # not(cloud_threshold) means that we're selecting all the intermediate values: anything between 
    # mask_clouds_lower and mask_clouds_upper.

    _cloud_threshold = (
        b7_landmasked_cloudmasked .< mask_clouds_lower .|| b7_landmasked_cloudmasked .> mask_clouds_upper
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
function clean_binary_floes(bw_img; min_opening_area=50)
    img_opened = area_opening(bw_img; min_area=min_opening_area) |> hbreak
    img_filled = branch(img_opened) |> bridge |> fill_holes
    diff_matrix = img_opened .!= img_filled
    return bw_img .|| diff_matrix
end

function segB_binarize(sharpened_image, brightened_image, cloudmask;
     gamma_factor=2.5, adjusted_ice_threshold=0.05, fill_range=(0, 1), alpha_level=0.5)
    # Weighted average between brightened image and sharpened grayscale
    adjusted_sharpened = (1 - alpha_level) .* sharpened_image .+ alpha_level .* brightened_image

    # Gamma correction and cloud masking
    adjust_histogram!(adjusted_sharpened, GammaCorrection(; gamma=gamma_factor))
    apply_cloudmask!(adjusted_sharpened, cloudmask) # Should this be happening here?

    # Thresholding and filling small holes (based on fill_range=(min, max))
    segB = adjusted_sharpened .<= adjusted_ice_threshold
    segb_filled = .!imfill(segB, fill_range)
    return segb_filled
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
    segmentation_F(
    brightened_image::Matrix{Gray{Float64}},
    segmentation_B_ice_intersect::BitMatrix,
    segmentation_B_watershed_intersect::BitMatrix,
    ice_labels::Union{Vector{Int64},BitMatrix},
    cloudmask::BitMatrix,
    landmask::BitMatrix;
    min_area_opening::Int64=20
)

Cleans up past segmentation images with morphological operations, and applies the results of prior watershed segmentation, returning the final cleaned image for tracking with ice floes segmented and isolated.

# Arguments
- `brightened_image`: Brightened grayscale image. Leads and interstitial ice should be dark. 
- `segmentation_B_ice_intersect`: binary mask output from `segmentation_b.jl`
- `segmentation_B_watershed_intersect`: ice pixels, output from `segmentation_b.jl`
- `ice_labels`: vector of pixel coordinates output from `find_ice_labels.jl`
- `cloudmask.jl`: bitmatrix cloudmask for region of interest
- `landmask.jl`: bitmatrix landmask for region of interest
- `min_area_opening`: threshold used for area opening; pixel groups greater than threshold are retained

"""
function segmentation_F( # rename
    brightened_image::Matrix{Gray{Float64}},
    segmentation_B_ice_intersect::BitMatrix,
    segmentation_B_watershed_intersect::BitMatrix,
    falsecolor_image,
    cloudmask::BitMatrix,
    landmask::BitMatrix;
    cluster_selection_algorithm,
    min_area_opening::Int64=20,
)::BitMatrix
    
# compare this workflow to the first instance with the k-means binarization.
    # Where were the cloudmask and landmasked used?    

    # segb_leads = 1 for leads/water, 0 for ice (and bright clouds)
    segb_leads = .!segmentation_B_watershed_intersect .* segmentation_B_ice_intersect
    segb_leads .= .!area_opening(segb_leads; min_area=min_area_opening, connectivity=2)
    
    # reconstruction by erosion is equivalent to dilation of the complements
    markers = dilate(brightened_image, strel_diamond((5, 5)))
    reconstructed_leads = complement.(mreconstruct(erode, markers, brightened_image)) 
    
    # multiply with the segb_leads, so if at least one of the two labels it as ice it stays ice
    reconstructed_leads .*= segb_leads 

    #### Update K-Means Segmentation ####

    leads_segmented =
        kmeans_binarization(reconstructed_leads, falsecolor_image;
            cluster_selection_algorithm=IceDetectionLopezAcosta2019()) .*
        .!segmentation_B_watershed_intersect

    return leads_segmented
end
"""

Series of morphological operations to produce final well-separated ice floes.

"""
function morphological_cleanup_final(
        binary_ice_img,
        cloudmask;
        fill_range=(0,1),
        min_area_opening=20,
        bothat_strel=strel_diamond((5,5)),
        opening_strel=centered(se_disk4()), # consider replacement with optimized strel 
    )
    
    img = deepcopy(binary_ice_img)
    img .= hbreak(img) 
    img .= branch(img)

    leads_filled = .!imfill(.!img, fill_range)
    leads_opened = area_opening(leads_filled; min_area=min_area_opening, connectivity=2)
    leads_opened .= branch(leads_opened)
    leads_bothat = bothat(leads_opened, bothat_strel) .> 0

    leads = (complement.(leads_bothat) .* leads_opened) .> 0
    area_opening!(leads, leads; min_area=min_area_opening, connectivity=2)
    # dmw: replace multiplication with apply_cloudmask
    leads_bothat_filled = (fill_holes(leads) .* .!cloudmask)
    
    floes = branch(leads_bothat_filled)
    floes_opened = opening(floes, opening_strel)
    mreconstruct!(dilate, floes_opened, floes, floes_opened)

    return floes_opened
end

"""
    _adjust_histogram(masked_view, nbins, rblocks, cblocks, clip)

Perform adaptive histogram equalization to a masked image. To be invoked within `imsharpen`.

# Arguments
- `masked_view`: input image in truecolor
See `imsharpen` for a description of the remaining arguments

"""
function _adjust_histogram(masked_view, nbins, rblocks, cblocks, clip)
    return adjust_histogram(
        masked_view,
        AdaptiveEqualization(;
            nbins=nbins,
            rblocks=rblocks,
            cblocks=cblocks,
            minval=minimum(masked_view),
            maxval=maximum(masked_view),
            clip=clip,
        ),
    )
end

"""
    imsharpen(truecolor_image, landmask_no_dilate, lambda, kappa, niters, nbins, rblocks, cblocks, clip, smoothing_param, intensity)

Sharpen `truecolor_image`.

# Arguments
- `truecolor_image`: input image in truecolor
- `landmask_no_dilate`: landmask for region of interest
- `lambda`: speed of diffusion (0–0.25)
- `kappa`: conduction coefficient for diffusion (25–100)
- `niters`: number of iterations of diffusion
- `nbins`: number of bins during histogram equalization
- `rblocks`: number of row blocks to divide input image during equalization
- `cblocks`: number of column blocks to divide input image during equalization
- `clip`: Thresholds for clipping histogram bins (0–1); values closer to one minimize contrast enhancement, values closer to zero maximize contrast enhancement
- `smoothing_param`: pixel radius for gaussian blurring (1–10)
- `intensity`: amount of sharpening to perform
"""
function imsharpen(
    truecolor_image::Matrix{RGB{Float64}},
    landmask_no_dilate::BitMatrix,
    lambda::Real=0.1,
    kappa::Real=0.1,
    niters::Int64=5,
    nbins::Int64=255,
    rblocks::Int64=10, # matlab default is 8 CP
    cblocks::Int64=10, # matlab default is 8 CP
    clip::Float64=0.86, # matlab default is 0.01 CP
    smoothing_param::Int64=10,
    intensity::Float64=2.0,
)::Matrix{Float64}
    input_image = apply_landmask(truecolor_image, landmask_no_dilate)

    pmd = PeronaMalikDiffusion(lambda, kappa, niters, "exponential")
    input_image .= nonlinear_diffusion(input_image, pmd)

    masked_view = Float64.(channelview(input_image))

    eq = [
        _adjust_histogram(@view(masked_view[i, :, :]), nbins, rblocks, cblocks, clip) for
        i in 1:3
    ]

    image_equalized = colorview(RGB, eq...)

    image_equalized_gray = Gray.(image_equalized)

    return unsharp_mask(image_equalized_gray, smoothing_param, intensity)
end

# TODO: Remove function, replace with direct use of landmask and colorview.
"""
    imsharpen_gray(imgsharpened, landmask)

Apply landmask and return Gray type image in colorview for normalization.

"""
function imsharpen_gray(
    imgsharpened::Matrix{Float64}, landmask::AbstractArray{Bool}
)::Matrix{Gray{Float64}}
    image_sharpened_landmasked = apply_landmask(imgsharpened, landmask)
    return colorview(Gray, image_sharpened_landmasked)
end


"""IceDetectionLopezAcosta2019

Application of the IceDetectionFirstNonZeroAlgorithm using two passes of 
the IceDetectionThresholdMODIS721 and one application of the IceDetectionBrightnessPeaksMODIS721.
"""
function IceDetectionLopezAcosta2019(;
    band_7_max::Float64=Float64(5 / 255),
    band_2_min::Float64=Float64(230 / 255),
    band_1_min::Float64=Float64(240 / 255),
    band_7_max_relaxed::Float64=Float64(10 / 255),
    band_1_min_relaxed::Float64=Float64(190 / 255),
    possible_ice_threshold::Float64=Float64(75 / 255),
)
    return IceDetectionFirstNonZeroAlgorithm([
        IceDetectionThresholdMODIS721(;
            band_7_max=band_7_max,
            band_2_min=band_2_min,
            band_1_min=band_1_min
        ),
        IceDetectionThresholdMODIS721(;
            band_7_max=band_7_max_relaxed,
            band_2_min=band_2_min,
            band_1_min=band_1_min_relaxed,
        ),
        IceDetectionBrightnessPeaksMODIS721(;
            band_7_max=band_7_max,
            possible_ice_threshold=possible_ice_threshold
        ),
    ])
end

end
