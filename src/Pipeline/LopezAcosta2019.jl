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
    stdmult,
    ⋅, # dot operator
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
    apply_cloudmask
import ..Segmentation:
    IceFloeSegmentationAlgorithm, 
    IceDetectionLopezAcosta2019,
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


@kwdef struct Segment <: IceFloeSegmentationAlgorithm
    landmask_structuring_element::AbstractMatrix{Bool} = make_landmask_se()
    cloud_mask_algorithm = LopezAcostaCloudMask(cloud_mask_thresholds...)
    diffusion_algorithm = PeronaMalikDiffusion(diffusion_parameters...)
    adapthisteq_params = (
        nbins=256,
        rblocks=10, # matlab default is 8 CP
        cblocks=10, # matlab default is 8 CP
        clip=0.86,  # matlab default is 0.01 CP
    )
    unsharp_mask_params = (smoothing_param=10, intensity=2)
    reconstruct_strel = strel_diamond((5, 5)) # structuring element used for enhancing foreground/background contrast
    kmeans_params = (k=4, maxiter=50, random_seed=45)
    ice_labels_algorithm = IceDetectionLopezAcosta2019(;ice_masks_params...)
end

# function Segment(; landmask_structuring_element=make_landmask_se())
#     return Segment(landmask_structuring_element)
# end

function (p::Segment)(
    truecolor::T,
    falsecolor::T,
    landmask::U;
    intermediate_results_callback::Union{Nothing,Function}=nothing,
) where {T<:AbstractMatrix{<:AbstractRGB},U<:AbstractMatrix}

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

    @info "Grayscale reconstruction of sharpened image"
    # input may need to be Gray.(sharpened_truecolor_image) in case the landmask dilation matters
    markers = complement.(dilate(sharpened_grayscale_image, p.reconstruct_strel))
    mask = complement.(sharpened_grayscale_image)
    # reconstruction of the complement: floes are dark, leads are bright
    reconstructed_grayscale = mreconstruct(dilate, markers, mask)
    apply_landmask!(reconstructed_grayscale, landmask_imgs.dilated)

    
    # Discriminate ice/water
    @info "Discriminating ice/water"
    ice_water_discrim = discriminate_ice_water(
        falsecolor_image, reconstructed_grayscale, copy(landmask_imgs.dilated), cloudmask
    )

    # 3. Segmentation
    @info "Segmenting floes part 1/3"
    # Components: 
    # - k-means binarization
    # - morphological cleanup
    # - apply cloudmask 
    # result labeled "segA" is a binarization image

    # Compute k-means clustering and detect ice masks
    # Important: Application of landmask before finding ice labels, to avoid at least some of the brightest landfast ice pixels.
    kmeans_result = kmeans_binarization(
        ice_water_discrim,
        apply_landmask(falsecolor_image, landmask_imgs.dilated);
        k=p.kmeans_params.k,
        maxiter=p.kmeans_params.maxiter,
        random_seed=p.kmeans_params.random_seed,
        ice_labels_algorithm=p.ice_labels_algorithm
        ) |> clean_binary_floes

    segA = apply_cloudmask(kmeans_result, cloudmask)

    # segmentation_B
    @info "Segmenting floes part 2/3"
    segB = segmentation_B(sharpened_grayscale_image, cloudmask, segA)

    # Process watershed in parallel using Folds
    @info "Building watersheds"
    watersheds_segB = [
        watershed_ice_floes(segB.not_ice_bit), watershed_ice_floes(segB.ice_intersect)
    ]
    watersheds_segB_product = watershed_product(watersheds_segB...)

    # segmentation_F
    # TODO: @hollandjg find out why segF is more dilated
    @info "Segmenting floes part 3/3"
    segF = segmentation_F(
        segB.not_ice,
        segB.ice_intersect,
        watersheds_segB_product,
        ice_mask,
        cloudmask,
        landmask_imgs.dilated,
    )

    @info "Labeling floes"
    labels = label_components(segF)

    # Return the original truecolor image, segmented
    segments = SegmentedImage(truecolor, labels)

    if !isnothing(intermediate_results_callback)
        segments_truecolor = SegmentedImage(truecolor, labels)
        segments_falsecolor = SegmentedImage(falsecolor, labels)
        intermediate_results_callback(;
            truecolor,
            falsecolor,
            landmask_dilated=landmask_imgs.dilated,
            landmask_non_dilated=landmask_imgs.non_dilated,
            cloudmask=cloudmask,
            ice_mask=ice_mask,
            sharpened_truecolor_image=sharpened_truecolor_image,
            sharpened_gray_truecolor_image=sharpened_grayscale_image,
            normalized_image=reconstructed_grayscale,
            ice_water_discrim=ice_water_discrim,
            segA=segA,
            segB=segB,
            watersheds_segB_product=watersheds_segB_product,
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
        falsecolor_image::Matrix{RGB{Float64}},
        normalized_image::Matrix{Gray{Float64}},
        landmask_bitmatrix::T,
        cloudmask_bitmatrix::T,
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
        nbins::Real=155
    )

Generates an image with ice floes apparent after filtering and combining previously processed versions
of falsecolor and truecolor images from the same region of interest. Returns an image ready for segmentation
to isolate floes.

# Arguments
- `falsecolor_image`: input image in false color reflectance
- `normalized_image`: normalized version of true color image
- `landmask_bitmatrix`: landmask for region of interest
- `cloudmask_bitmatrix`: cloudmask for region of interest
- `floes_threshold`: heuristic applied to original false color image
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
    falsecolor_image::Matrix{RGB{Float64}},
    normalized_image::Matrix{Gray{Float64}},
    landmask_bitmatrix::T,
    cloudmask_bitmatrix::T,
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
    nbins::Real=155,
)::AbstractMatrix where {T<:AbstractArray{Bool}}
    clouds_channel = create_clouds_channel(cloudmask_bitmatrix, falsecolor_image)
    falsecolor_image_band7 = @view(channelview(falsecolor_image)[1, :, :])

    # first define all of the image variations
    image_clouds = apply_landmask(clouds_channel, landmask_bitmatrix) # output during cloudmask apply, landmasked
    image_cloudless = apply_landmask(falsecolor_image_band7, landmask_bitmatrix) # channel 1 (band 7) from source falsecolor image, landmasked
    image_floes = apply_landmask(falsecolor_image, landmask_bitmatrix) # source false color reflectance, landmasked
    image_floes_view = channelview(image_floes)

    floes_band_2 = @view(image_floes_view[2, :, :])
    floes_band_1 = @view(image_floes_view[3, :, :])

    # keep pixels greater than intensity 100 in bands 2 and 1
    floes_band_2_keep = floes_band_2[floes_band_2 .> floes_threshold]
    floes_band_1_keep = floes_band_1[floes_band_1 .> floes_threshold]

    _, floes_bin_counts = build_histogram(floes_band_2_keep, nbins)
    _, vals = findmaxima(floes_bin_counts)

    differ = vals / (maximum(vals))
    proportional_intensity = sum(differ .> differ_threshold) / length(differ) # finds the proportional intensity of the peaks in the histogram

    # compute kurtosis, skewness, and standard deviation to use in threshold filtering
    kurt_band_2 = kurtosis(floes_band_2_keep)
    skew_band_2 = skewness(floes_band_2_keep)
    kurt_band_1 = kurtosis(floes_band_1_keep)
    standard_dev = stdmult(⋅, normalized_image)

    # find the ratio of clouds in the image to use in threshold filtering
    _, clouds_bin_counts = build_histogram(image_clouds .> 0)
    total_clouds = sum(clouds_bin_counts[51:end])
    total_all = sum(clouds_bin_counts)
    clouds_ratio = total_clouds / total_all

    threshold_50_check = _check_threshold_50(
        kurt_band_1,
        kurt_band_2,
        kurt_thresh_lower,
        kurt_thresh_upper,
        skew_band_2,
        skew_thresh,
        proportional_intensity,
    )

    threshold_130_check = _check_threshold_130(
        clouds_ratio,
        clouds_ratio_threshold,
        standard_dev,
        st_dev_thresh_lower,
        st_dev_thresh_upper,
    )

    if threshold_50_check
        THRESH = 50 / 255
    elseif threshold_130_check
        THRESH = 130 / 255
    else
        THRESH = 80 / 255 #intensity value of 80
    end

    normalized_image_copy = copy(normalized_image)
    normalized_image_copy[normalized_image_copy .> THRESH] .= 0
    @. normalized_image_copy = normalized_image - (normalized_image_copy * 3)

    lm = deepcopy(landmask_bitmatrix)
    @. lm = (image_clouds < mask_clouds_lower || image_clouds > mask_clouds_upper)
    @. image_cloudless = image_cloudless * !lm
    @. normalized_image_copy = clamp01nan(normalized_image_copy - (image_cloudless * 3))

    return normalized_image_copy
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
    struct_elem=strel_diamond((3, 3));
    fill_range::Tuple=(0, 1),
    isolation_threshold::Float64=0.4,
    alpha_level::Float64=0.5,
    gamma_factor::Float64=2.5,
    adjusted_ice_threshold::Float64=0.05,
)

    ## Process sharpened image
    not_ice_mask = deepcopy(sharpened_image)
    not_ice_mask[not_ice_mask .< isolation_threshold] .= 0
    not_ice_bit = not_ice_mask .* 0.3
    not_ice_mask .= not_ice_bit .+ sharpened_image
    adjusted_sharpened = (
        (1 - alpha_level) .* sharpened_image .+ alpha_level .* not_ice_mask
    )
    gamma_adjusted_sharpened = adjust_histogram(
        adjusted_sharpened, GammaCorrection(; gamma=gamma_factor)
    )
    gamma_adjusted_sharpened_cloudmasked = apply_cloudmask(
        gamma_adjusted_sharpened, cloudmask
    )
    segb_filled =
        .!imfill(
            gamma_adjusted_sharpened_cloudmasked .<= adjusted_ice_threshold, fill_range
        )

    ## Process ice mask
    segb_ice = closing(segmented_a_ice_mask, struct_elem) .* segb_filled

    ice_intersect = (segb_filled .* segb_ice)

    return (;
        :not_ice => map(clamp01nan, not_ice_mask)::Matrix{Gray{Float64}},
        :not_ice_bit => (not_ice_bit .> 0)::BitMatrix,
        :ice_intersect => ice_intersect::BitMatrix,
    )
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
    segmentation_F(
    segmentation_B_not_ice_mask::Matrix{Gray{Float64}},
    segmentation_B_ice_intersect::BitMatrix,
    segmentation_B_watershed_intersect::BitMatrix,
    ice_labels::Union{Vector{Int64},BitMatrix},
    cloudmask::BitMatrix,
    landmask::BitMatrix;
    min_area_opening::Int64=20
)

Cleans up past segmentation images with morphological operations, and applies the results of prior watershed segmentation, returning the final cleaned image for tracking with ice floes segmented and isolated.

# Arguments
- `segmentation_B_not_ice_mask`: gray image output from `segmentation_b.jl`
- `segmentation_B_ice_intersect`: binary mask output from `segmentation_b.jl`
- `segmentation_B_watershed_intersect`: ice pixels, output from `segmentation_b.jl`
- `ice_labels`: vector of pixel coordinates output from `find_ice_labels.jl`
- `cloudmask.jl`: bitmatrix cloudmask for region of interest
- `landmask.jl`: bitmatrix landmask for region of interest
- `min_area_opening`: threshold used for area opening; pixel groups greater than threshold are retained

"""
function segmentation_F(
    segmentation_B_not_ice_mask::Matrix{Gray{Float64}},
    segmentation_B_ice_intersect::BitMatrix,
    segmentation_B_watershed_intersect::BitMatrix,
    ice_labels::Union{Vector{Int64},BitMatrix,AbstractArray{<:Gray}},
    cloudmask::BitMatrix,
    landmask::BitMatrix;
    min_area_opening::Int64=20,
)::BitMatrix
    apply_landmask!(segmentation_B_not_ice_mask, landmask)

    ice_leads = .!segmentation_B_watershed_intersect .* segmentation_B_ice_intersect

    ice_leads .= .!area_opening(ice_leads; min_area=min_area_opening, connectivity=2)

    not_ice = dilate(segmentation_B_not_ice_mask, strel_diamond((5, 5)))

    mreconstruct!(
        dilate, not_ice, complement.(not_ice), complement.(segmentation_B_not_ice_mask)
    )

    reconstructed_leads = (not_ice .* ice_leads) .+ (60 / 255)

    leads_segmented =
        kmeans_segmentation(reconstructed_leads, ice_labels) .*
        .!segmentation_B_watershed_intersect
    @info("Done with k-means segmentation")
    leads_segmented_broken = hbreak(leads_segmented)

    leads_branched = branch(leads_segmented_broken)

    leads_filled = .!imfill(.!leads_branched, 0:1)

    leads_opened = branch(
        area_opening(leads_filled; min_area=min_area_opening, connectivity=2)
    )

    leads_bothat = bothat(leads_opened, strel_diamond((5, 5))) .> 0.499

    leads = convert(BitMatrix, (complement.(leads_bothat) .* leads_opened))

    area_opening!(leads, leads; min_area=min_area_opening, connectivity=2)

    # dmw: replace multiplication with apply_cloudmask
    leads_bothat_filled = (fill_holes(leads) .* .!cloudmask)
    # leads_bothat_filled = apply_cloudmask(fill_holes(leads), cloudmask)
    floes = branch(leads_bothat_filled)

    floes_opened = opening(floes, centered(se_disk4()))

    mreconstruct!(dilate, floes_opened, floes, floes_opened)

    return floes_opened
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

    # Diffusion algorithm can be an input to the overall algorithm
    pmd = PeronaMalikDiffusion(lambda, kappa, niters, "exponential")
    input_image .= nonlinear_diffusion(input_image, pmd)

    # Channelview -> _adjusthistogram -> colorview as a function.
    masked_view = Float64.(channelview(input_image))

    eq = [
        _adjust_histogram(@view(masked_view[i, :, :]), nbins, rblocks, cblocks, clip) for
        i in 1:3
    ]

    image_equalized = colorview(RGB, eq...)

    image_equalized_gray = Gray.(image_equalized)

    return unsharp_mask(image_equalized_gray, smoothing_param, intensity)
end

end