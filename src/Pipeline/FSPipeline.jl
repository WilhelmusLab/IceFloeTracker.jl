module FSPipeline
"""
Update of the LopezAcosta2019 algorithm for the floe trajectory dataset
"""

using Images
import Peaks: findmaxima
import StatsBase: kurtosis, skewness, mean, std

import ..Filtering: nonlinear_diffusion, PeronaMalikDiffusion, unsharp_mask, channelwise_adapthisteq
import ..Morphology: hbreak, hbreak!, branch, bridge, fill_holes, se_disk4
import ..Preprocessing:
    create_landmask,
    create_cloudmask,
    apply_landmask,
    apply_landmask!,
    apply_cloudmask,
    apply_cloudmask!,
    LopezAcostaCloudMask,
    Watkins2025CloudMask
import ..ImageUtils: get_tiles
import ..Segmentation: 
    IceFloeSegmentationAlgorithm, 
    find_ice_mask, 
    kmeans_binarization,
    IceDetectionFirstNonZeroAlgorithm,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionThresholdMODIS721,
    stitch_clusters

"""
Sample input parameters expected by the main function
"""
# cloud_mask_thresholds = (
#     prelim_threshold=110.0 / 255.0,
#     band_7_threshold=200.0 / 255.0,
#     band_2_threshold=190.0 / 255.0,
#     ratio_lower=0.0,
#     ratio_offset=0.0,
#     ratio_upper=0.75,
# )

diffusion_parameters = (lambda=0.1, kappa=0.1, niters=5, g="exponential")

@kwdef struct Segment <: IceFloeSegmentationAlgorithm
    coastal_buffer_structuring_element::AbstractMatrix{Bool} = strel_box((51,51))
    cloud_mask_algorithm = Watkins2025CloudMask()
    diffusion_algorithm = PeronaMalikDiffusion(diffusion_parameters...)
    adapthisteq_params = (
        nbins=256,
        rblocks=16, # TODO: Make responsive to tile size
        cblocks=8, 
        clip=0.9,  # matlab default is 0.01 CP, which should be the same as clip=0.99
    )
    tile_size_pixels=1000
    min_tile_ice_pixel_count=500
    unsharp_mask_params = (smoothing_param=10, intensity=0.5)
    kmeans_params = (k=4, maxiter=50, random_seed=45)
    cluster_selection_algorithm = IceDetectionBrightnessPeaksMODIS721(
        band_7_max=0.1,
        possible_ice_threshold=0.3,
        join_method="union",
        minimum_prominence=0.01)
    grayscale_floe_threshold = 0.4 # Brighten pixels brighter than this threshold 
    adaptive_binarization_settings = (
        minimum_window_size=400,
        threshold_percentage=0,
        minimum_brightness=100/255) # Not used in this workflow yet
    watershed_strel = se_disk(5) # 
end

function (p::Segment)(
    truecolor::T,
    falsecolor::T,
    land_mask::U;
    intermediate_results_callback::Union{Nothing,Function}=nothing,
) where {T<:AbstractMatrix{<:AbstractRGB},U<:AbstractMatrix}
    @info "building landmask and coastal buffer mask"
    _landmask_temp = create_landmask(
        float64.(land_mask), p.coastal_buffer_structuring_element
    )
    land_mask = _landmask_temp.non_dilated
    coastal_buffer = _landmask_temp.dilated
    return p(
        truecolor,
        falsecolor,
        land_mask,
        coastal_buffer;
        intermediate_results_callback=intermediate_results_callback,
    )
end

function (p::Segment)(
    truecolor::T,
    falsecolor::T,
    land_mask::U,
    coastal_buffer::U;
    intermediate_results_callback::Union{Nothing,Function}=nothing,
) where {T<:AbstractMatrix{<:AbstractRGB},U<:BitMatrix}

    # Move these conversions down through the function as each step gets support for 
    # the full range of image formats
    truecolor_image = float64.(RGB.(truecolor))
    falsecolor_image = float64.(RGB.(falsecolor))

    n, m = size(truecolor_image)
    tile_size_pixels = p.tile_size_pixels
    tile_size_pixels > maximum([n, m]) && begin
        @warn "Tile size too large, defaulting to image size"
            tile_size_pixels = minimum([n, m])
    end
    tiles = get_tiles(truecolor_image, tile_size_pixels)
 
    @info "Building masks"
    # TODO: Make sure tests aren't over-sensitive to roundoff errors for Float32 vs Float64
    cloud_mask = create_cloudmask(falsecolor_image, p.cloud_mask_algorithm)

    # 2. Intermediate images
    apply_landmask!(falsecolor_image, coastal_buffer)
    apply_cloudmask!(falsecolor_image, cloud_mask)

    prelim_ice_mask = (red.(truecolor_image) .> 0.3) .&& .!cloud_mask .&& .!coastal_buffer
    filtered_tiles = filter(
        t -> sum(prelim_ice_mask[t...]) > p.min_tile_ice_pixel_count, # check if the min is larger than the tile size and adjust
        tiles);

    @info "Preprocessing truecolor image"
    apply_landmask!(truecolor_image, land_mask);
    preproc_img = float64.(deepcopy(truecolor_image));

    # TODO: Options to test
    # Order of operations permutations: 
    #    a. eq -> sharpen -> diffuse
    #    b. diffuse -> eq -> sharpen (current Lopez Acosta 2019)
    #    c. sharpen -> diffuse -> eq
    #    d. eq -> diffuse -> sharpen
    #    e. diffuse -> sharpen -> eq
    #    f. sharpen -> diffuse -> eq
    #
    # Color image processing vs grayscale image processing
    #    a. Grayscale immediately
    #    b. Grayscale before sharpening (as is currently)
    #    c. Update sharpen to be channel-wise, then cast to grayscale afterward
    for tile in filtered_tiles
        preproc_img[tile...] .= nonlinear_diffusion(
            preproc_img[tile...], p.diffusion_algorithm
        )
    end
    
    preproc_img .= channelwise_adapthisteq(preproc_img;
        nbins=p.adapthisteq_params.nbins,
        rblocks=p.adapthisteq_params.rblocks,
        cblocks=p.adapthisteq_params.cblocks,
        clip=p.adapthisteq_params.clip
    );

    preproc_gray = Gray.(preproc_img)
    for tile in filtered_tiles
        preproc_gray[tile...] .= unsharp_mask(
        preproc_gray[tile...],
        p.unsharp_mask_params.smoothing_param,
        p.unsharp_mask_params.intensity)
    end

    apply_landmask!(preproc_gray, coastal_buffer)
    apply_cloudmask!(preproc_gray, cloud_mask);

    
    @info "Enhancing grayscale image"
    # This step is a grayscale morphology operation. Reconstruction by dilation of the image complement
    # followed by thresholding.
    ice_water_discrim = zeros(size(preproc_gray))
    for tile in filtered_tiles
        ice_water_discrim[tile...] .= discriminate_ice_water(
            preproc_gray[tile...], falsecolor_image[tile...], coastal_buffer[tile...], cloud_mask[tile...]
    );
    end
 
    # 3. Segmentation
    @info "Segmenting floes part 1/3"
    kmeans_result = kmeans_binarization(
            Gray.(ice_water_discrim),
            falsecolor_image,
            filtered_tiles;
            k=p.kmeans_params.k,
            maxiter=p.kmeans_params.maxiter,
            random_seed=p.kmeans_params.random_seed,
            cluster_selection_algorithm=p.cluster_selection_algorithm
            ) |> clean_binary_floes
    
    # check: are there any regions that are nonzero under the cloudmask, since it was applied in discriminate ice water?
    apply_cloudmask!(kmeans_result, cloud_mask)

    # The clean binary floes method has an aggressive fill_holes algorithm. Potentially merging with the
    # ice brightness threshold can prevent some of the interstitial water areas from being filled.

    # segmentation_B
    @info "Segmenting floes part 2/3"
    # TODO: Bring over the updated version I've been developing for the streamline algorithm
    # grayscale enhnacement using gamma correction
    # threshold binarization
    # morphological cleanup

    segB = segmentation_B(preproc_gray, cloud_mask, kmeans_result)

    # Process watershed in parallel using Folds
    @info "Computing watershed boundaries"
    w_merged = watershed_transform(
        segB.ice_intersect,
        segB.not_ice,
        filtered_tiles;
        strel=s.watershed_strel,
        dist_threshold=4
    )
    
    @info "tiled watershed on the not ice bit, whatever that one is"
    w_other = watershed_transform(
        segB.not_ice_bit,
        segB.not_ice,
        filtered_tiles;
        strel=s.watershed_strel, # add to inputs for calibration
        dist_threshold=4 # add to inputs for calibration
    )

    watersheds_segB_product = falses(size(truecolor_image))
    for tile in filtered_tiles
        watersheds_segB_product[tile...] .= (isboundary(labels_map(w_merged)[tile...]) .* isboundary(labels_map(w_other)[tile...])) .> 0
    end

    
    # segmentation_F
    # TODO: Split the refined k-means workflow from the cleanup
    # TODO: Make the 
    @info "Segmenting floes part 3/3"
    segF = segmentation_F(
        segB.not_ice,
        segB.ice_intersect,
        watersheds_segB_product,
        falsecolor_image,
        cloud_mask,
        coastal_buffer;
        cluster_selection_algorithm=p.cluster_selection_algorithm
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
            coastal_buffer,
            cloud_mask,
            prelim_ice_mask,
            sharpened_grayscale_image=preproc_gray,
            ice_water_discrim=ice_water_discrim,
            segA=kmeans_result,
            segB=segB,
            watersheds_segB_product=watersheds_segB_product,
            segF=segF,
            labels=labels,
            segment_mean_truecolor=map( # TODO Add "view_seg" code snippet
                i -> segment_mean(segments_truecolor, i), labels_map(segments_truecolor)
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
    differ_threshold::Float64=0.6
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
    # clouds lower, or brighter than clouds lower.
    # So the b7_landmasked getting multiplied by not(cloud_threshold) means that we're selecting all
    # the intermediate values: anything between mask_clouds_lower and mask_clouds_upper.
    # This is essentially speckle noise, and has minimal effect on the results.

    # _cloud_threshold = (
    #     b7_landmasked_cloudmasked .< mask_clouds_lower .|| b7_landmasked_cloudmasked .> mask_clouds_upper
    # )

    # reusing image_cloudless - used to be band7_masked
    # I think the names were backwards: image_cloudless had not been cloudmasked, and image_clouds had been.
    # So the question is if they're swapped in 305 and 311 also.
    # @. b7_landmasked = b7_landmasked * !_cloud_threshold

    # Check to see if selecting indices to set to 0 would be equivalent
    # Also check if rescaling intensity would be safer
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

"""
 segB_binarize(sharpened_image, brightened_image, cloudmask;
     gamma_factor=2.5, adjusted_ice_threshold=0.05, fill_range=(0, 1), alpha_level=0.5)

Binarize the sharpened image by selective brightening, gamma correction, threshold application, and 
clean up with image hole filling.

"""
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


# """
#     segmented_ice_cloudmasking(gray_image, cloudmask, ice_labels;)

# Apply cloudmask to a bitmatrix of segmented ice after kmeans clustering. Returns a bitmatrix with open water/clouds = 0, ice = 1).

# # Arguments

# - `gray_image`: output image from `ice-water-discrimination.jl` or gray ice floe leads image in `segmentation_f.jl`
# - `cloudmask`: bitmatrix cloudmask for region of interest
# - `ice_labels`: vector if pixel coordinates output from `find_ice_labels.jl`

# """
# function segmented_ice_cloudmasking(
#     gray_image, 
#     falsecolor_image,
#     cloudmask::BitMatrix,
# )::BitMatrix
#     segmented_ice = kmeans_binarization(
#         gray_image, falsecolor_image; 
#         cluster_selection_algorithm=IceDetectionLopezAcosta2019())
#     segmented_ice_cloudmasked = deepcopy(segmented_ice)
#     segmented_ice_cloudmasked[cloudmask] .= 0
#     return segmented_ice_cloudmasked
# end

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

""" # TODO: Split into the segB binarize and segB preproc functions
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

##### Watershed transformation #######
"""
   watershed_transform(binary_img, img; strel=strel_diamond((3,3)), dist_threshold=4)
   watershed_transform(binary_img, img, tiles; strel=strel_diamond((3,3)), dist_threshold=4)

    Carry out a watershed transform of `binary_img` by computing the distance transform,
    selecting marker regions with distance to background greater than `dist_threshold`, then eroding
    the marker regions using `strel`. Following the watershed transform, the background of `binary_img`
    is imposed and a SegmentedImage is generated using `img`. Optionally, supply a TiledIterator and 
    carry out the transform only on the listed tiles.

""" # TODO: Make the "marker_selection_function" an input, with the idea that it's a function of the distance transform
function watershed_transform(
    binary_floes,
    img;
    strel=strel_diamond((5,5)),
    dist_threshold=4
)
    bw = .!binary_floes
    dist = .- distance_transform(feature_transform(bw))
    markers = erode(dist .< -dist_threshold, strel) |> label_components
    labels = labels_map(watershed(dist, markers))
    labels[bw] .= 0

    return SegmentedImage(img, labels)
end

function watershed_transform(
    binary_floes,
    img,
    tiles;
    strel=strel_diamond((5,5)),
    dist_threshold=4,
    min_overlap=10,
    grayscale_threshold=0.1
)
    labels = zeros(Int64, size(binary_floes))
    for tile in tiles
        wseg = watershed_transform(binary_floes[tile...], img[tile...]; strel=strel, dist_threshold=dist_threshold)
        labels[tile...] = labels_map(wseg)
    end
    wseg = SegmentedImage(img, labels)
    stitched_labels = stitch_clusters(wseg, tiles, min_overlap, grayscale_threshold) # TODO: Add method to allow labels map for stitch clusters
    stitched_labels[.!binary_floes] .= 0 # TODO: Check to see if the background needs to be re-masked

    return SegmentedImage(img, stitched_labels)
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

""" # TODO: Split into the grayscale and floe splitting components
function segmentation_F(
    segmentation_B_not_ice_mask::Matrix{Gray{Float64}},
    segmentation_B_ice_intersect::BitMatrix,
    segmentation_B_watershed_intersect::BitMatrix,
    falsecolor_image,
    cloudmask::BitMatrix,
    landmask::BitMatrix;
    min_area_opening::Int64=20,
    cluster_selection_algorithm
)::BitMatrix
    apply_landmask!(segmentation_B_not_ice_mask, landmask)

    ice_leads = .!segmentation_B_watershed_intersect .* segmentation_B_ice_intersect

    ice_leads .= .!area_opening(ice_leads; min_area=min_area_opening, connectivity=2)

    not_ice = dilate(segmentation_B_not_ice_mask, strel_diamond((5, 5)))

    mreconstruct!(
        dilate, not_ice, complement.(not_ice), complement.(segmentation_B_not_ice_mask)
    )

    reconstructed_leads = (not_ice .* ice_leads) .+ (60 / 255)

    #### Update K-Means Segmentation ####
    @info("Updated k-means segmentation")
    leads_segmented =
        kmeans_binarization(reconstructed_leads, falsecolor_image;
            cluster_selection_algorithm=cluster_selection_algorithm) .*
        .!segmentation_B_watershed_intersect
    
    @info "Splitting floes"
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

    floes_opened = opening(floes, se_disk(4))

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
    se_disk(r)

    Generate an approximately circular structuring element with radius r. For small r, this will be somewhat diamond-shaped.

""" # #TODO add simple example to docs, add to special strels, and in future, optimize the extreme filter for this shape
function se_disk(r)
    se = [sum(abs.(c.I .- (r + 1)) .^ 2) for c in CartesianIndices((2*r + 1, 2*r + 1))]
    return sqrt.(se) .<= r
end

end
