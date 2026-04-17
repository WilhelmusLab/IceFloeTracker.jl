module FSPipeline
"""
Update of the LopezAcosta2019 algorithm for the Fram Strait floe trajectory dataset. Key updates:
    * Separating segmentation algorithm components for clearer functions
    * Implementation of the tiling workflow
    * Simpler watershed transformation steps
    * Exposure of algorithm settings to the user
"""

using Images
import Peaks: findmaxima
import StatsBase: kurtosis, skewness, mean, std

import ..Filtering: 
    nonlinear_diffusion, 
    PeronaMalikDiffusion, 
    unsharp_mask, 
    ContrastLimitedAdaptiveHistogramEqualization

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
import ..ImageUtils: get_tiles, imbrighten
import ..Segmentation: 
    IceFloeSegmentationAlgorithm, 
    find_ice_mask, 
    kmeans_binarization,
    tiled_adaptive_binarization,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionBrightnessPeaksMODIS134,
    stitch_clusters,
    view_seg

abstract type IceFloePreprocessingAlgorithm end


"""
   Preprocess(
        coastal_buffer_structuring_element = strel_box((51,51))
        cloud_mask_algorithm = Watkins2025CloudMask()
        diffusion_algorithm = PeronaMalikDiffusion(lambda=0.1, kappa=0.1, niters=5, g="exponential")
        adapthisteq_params = (nbins=256, rblocks=8, cblocks=8, clip=0.95)
        unsharp_mask_params = (radius=50, amount=0.2, threshold=0.01)
    )
    Preprocess()(img, cloudmask, landmask)

    Converts input image to grayscale, then preprocesses by appling nonlinear diffusion, 
    adaptive histogram equalization, and unsharp masking. Diffusion and unsharp masking are applied 
    to each tile, while the adaptive histogram equalization is divided according to the parameter specifications.

"""
@kwdef struct Preprocess <: IceFloePreprocessingAlgorithm
    diffusion_algorithm = PeronaMalikDiffusion(λ=0.1, K=0.1, niters=5, g="exponential")
    adapthisteq_params = (nbins=256, rblocks=8, cblocks=8, clip=0.99) # rblocks/cblocks not used yet -- add with CLAHE.jl
    unsharp_mask_params = (radius=50, amount=0.2, threshold=0.01)
end

function (p::Preprocess)(
    truecolor_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}}, 
    landmask,
    cloud_mask,
    tiles
)
    # Cast to grayscale first to save compute time
    proc_img = Gray.(truecolor_image)
    apply_landmask!(proc_img, landmask .|| cloud_mask)

    # Diffusion and sharpening
    proc_img .= nonlinear_diffusion(proc_img, tiles, p.diffusion_algorithm)
    
    adjust_histogram!(proc_img,
        ContrastLimitedAdaptiveHistogramEqualization(
            nbins=p.adapthisteq_params.nbins,
            rblocks=p.adapthisteq_params.rblocks,
            cblocks=p.adapthisteq_params.cblocks,
            clip=p.adapthisteq_params.clip)
    )

    for tile in tiles
        proc_img[tile...] .= unsharp_mask(proc_img[tile...],
                p.unsharp_mask_params.radius,
                p.unsharp_mask_params.amount,
                p.unsharp_mask_params.threshold)
    end

            
    apply_landmask!(proc_img, landmask .|| cloud_mask)
    
    return proc_img
end

"""
    FSPipeline.Segment()

Segmentation routine for identifying moderate to large floes in the Fram Strait. 
The image preprocessing is supplied as an function in the functor setup.


# Parameters
- `coastal_buffer_structuring_element::AbstractMatrix{Bool} = strel_box((51,51))`: Structuring element for the `create_landmask` function
- `cloud_mask_algorithm = LopezAcostaCloudMask()`: Cloud mask algorithm
- `preprocessing_algorithm = Preprocess()`: Function to sharpen and equalize the truecolor image
- `tile_size_pixels=1000`: Nominal tile size in pixels
- `min_tile_ice_pixel_count=300`: Smallest number of likely sea ice floes in tile
- `min_floe_size=100`: Smallest floe to retain in image
- `preliminary_ice_mask = IceDetectionBrightnessPeaksMODIS134(band_7_max=0.1, possible_ice_threshold=0.3)`: Function to use to identify likely ice pixels for filtering.
- `kmeans_params = (k=4, maxiter=50, random_seed=45)`: Parameters for `kmeans_binarization`
- `cluster_selection_algorithm = IceDetectionBrightnessPeaksMODIS721(
    band_7_max=0.1,
    possible_ice_threshold=0.3,
    join_method="union",
    minimum_prominence=0.01)`: Function to use to select a k-means cluster in the `kmeans_binarization` workflow
- `brightening_factor = 0.3`
- `adaptive_binarization_settings = (
    minimum_window_size=400,
    threshold_percentage=0,
    minimum_brightness=100/255)`: Settings for the adaptive threshold binarization
- `watershed_strel = se_disk(5)`
- `floe_splitting_settings = (max_fill_area=1, min_area_opening=20, opening_strel=se_disk(2))`
"""
@kwdef struct Segment <: IceFloeSegmentationAlgorithm # Tried making this an IceFloeSegmentationAlgorithm but julia complains about ambiguity when I do so. Need to fix this to be able to use the run and validate function.
    coastal_buffer_structuring_element::AbstractMatrix{Bool} = strel_box((51,51))
    cloud_mask_algorithm = LopezAcostaCloudMask()
    preprocessing_algorithm = Preprocess()
    tile_size_pixels=1000
    min_tile_ice_pixel_count=300
    min_floe_size=100
    kmeans_params = (k=4, maxiter=50, random_seed=45)
    preliminary_ice_mask = IceDetectionBrightnessPeaksMODIS134(band_1_min=0.3)
    cluster_selection_algorithm = IceDetectionBrightnessPeaksMODIS721(
        band_7_max=0.1,
        possible_ice_threshold=0.3,
        join_method="union",
        minimum_prominence=0.01)
    brightening_factor = 0.3 # Fraction to brighten floes by (e.g., 0 = no brightening, 1 = doubling brightness)
    adaptive_binarization_settings = (
        minimum_window_size=400,
        threshold_percentage=0,
        minimum_brightness=100/255)
    watershed_strel = se_disk(5) # 
    floe_splitting_settings = (max_fill_area=1, min_area_opening=20, opening_strel=se_disk(2))
end 

function (p::Segment)(
    truecolor::T₁,
    falsecolor::T₂,
    landmask::T₃,
    coastal_buffer_mask::T₄;
    intermediate_results_callback::Union{Nothing,Function}=nothing,
) where {
    T₁<:AbstractMatrix{<:Union{AbstractRGB,TransparentRGB}},
    T₂<:AbstractMatrix{<:Union{AbstractRGB,TransparentRGB}},
    T₃<:AbstractMatrix{<:Union{Bool,Gray{Bool}}},
    T₄<:AbstractMatrix{<:Union{Bool,Gray{Bool}}},
}
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
    # TODO: Option to supply number of row and cols 
    tiles = get_tiles(truecolor_image, tile_size_pixels)
 
    @info "Building masks"
    # TODO: Make sure tests aren't over-sensitive to roundoff errors for Float32 vs Float64
    cloud_mask = create_cloudmask(falsecolor_image, p.cloud_mask_algorithm)
    landmask = landmask .> 0 # make sure it's a bitmatrix
    # 2. Intermediate images - using coastal buffer on the FC image
    apply_landmask!(truecolor_image, landmask .|| cloud_mask)
    apply_landmask!(falsecolor_image, coastal_buffer_mask .|| cloud_mask)

    prelim_ice_mask = p.preliminary_ice_mask(truecolor_image, tiles)
    filtered_tiles = filter(
        t -> sum(prelim_ice_mask[t...]) > p.min_tile_ice_pixel_count, tiles);

    @info "Preprocessing truecolor image"
    preproc_gray = float64.(p.preprocessing_algorithm(
        truecolor_image, landmask, cloud_mask, filtered_tiles));
    
    # @info "Enhancing grayscale image"
    # # This step is a grayscale morphology operation. Reconstruction by dilation of the image complement
    # # followed by thresholding.
    # ice_water_discrim = zeros(size(preproc_gray))
    # for tile in filtered_tiles
    #     ice_water_discrim[tile...] .= discriminate_ice_water(
    #         preproc_gray[tile...], falsecolor_image[tile...], coastal_buffer_mask[tile...], cloud_mask[tile...]
    # );
    # end
 
    # 3. Segmentation
    @info "Segmenting floes part 1/3"
    # Compare to Segmentation A from LopezAcosta2019
    kmeans_result = kmeans_binarization(
            Gray.(preproc_gray),
            falsecolor_image, 
            filtered_tiles;
            k=p.kmeans_params.k,
            maxiter=p.kmeans_params.maxiter,
            random_seed=p.kmeans_params.random_seed,
            cluster_selection_algorithm=p.cluster_selection_algorithm
            )
    # Simpler cleanup
    # kmeans_result = opening(kmeans_result, se_disk(2))
    # kmeans_result .= imfill(kmeans_result, (0, p.min_floe_size))

    kmeans_result .= clean_binary_floes2(kmeans_result, prelim_ice_mask, cloud_mask)

    
    # check: are there any regions that are nonzero under the cloudmask, since it was applied in discriminate ice water?
    apply_landmask!(kmeans_result, landmask)
    apply_cloudmask!(kmeans_result, cloud_mask)

    # The clean binary floes method has an aggressive fill_holes algorithm. Potentially merging with the
    # ice brightness threshold can prevent some of the interstitial water areas from being filled.
    
    @info "Segmenting floes part 2/3"
    # Compare to Segmentation B from Lopez-Acosta 2019
    # grayscale enhnacement using gamma correction - does it ever help? Or is it always just super bright and all ice?

    # segB grayscale enhancement
    adaptive_binarized =  tiled_adaptive_binarization(preproc_gray, filtered_tiles; p.adaptive_binarization_settings...) .> 0 
    brightened_gray = imbrighten(preproc_gray, adaptive_binarized, 1 + p.brightening_factor)
    
    # # TODO: give more descriptive name
    # gamma_binarized = segB_binarize(preproc_gray, brightened_gray, cloud_mask; 
    #                                 gamma_factor=2.5, 
    #                                 adjusted_ice_threshold=0.05,
    #                                 fill_range=(0, 1),
    #                                 alpha_level=0.5)

    ice_intersect = closing(kmeans_result, strel_diamond((3, 3))) # .* gamma_binarized
    
    # Process watershed in parallel using Folds
    @info "Computing watershed boundaries"
    w_merged = watershed_transform(
        ice_intersect,
        brightened_gray,
        filtered_tiles;
        strel=p.watershed_strel,
        dist_threshold=4
    )
    
    w_other = watershed_transform(
        adaptive_binarized,
        brightened_gray,
        filtered_tiles;
        strel=p.watershed_strel,
        dist_threshold=4 # add to inputs for calibration
    )

    watersheds_product = falses(size(truecolor_image))
    for tile in filtered_tiles
        watersheds_product[tile...] .= (
            isboundary(labels_map(w_merged)[tile...]) .* isboundary(labels_map(w_other)[tile...])
            ) .> 0
    end


    # segmentation_F
    # TODO: Split the refined k-means workflow from the cleanup
    @info "Segmenting floes part 3/3"
    morphed_grayscale = reconstruct_and_mask(
        brightened_gray,
        watersheds_product,
        ice_intersect,
        landmask
    ) # allows min area and strel to be added as kwd inputs

    # kmeans binarization, again
    segF_binarized = kmeans_binarization(
        morphed_grayscale, falsecolor_image, filtered_tiles; # number of clusters = ?
        cluster_selection_algorithm=p.cluster_selection_algorithm
        ) .* .! watersheds_product
    
    @info "Splitting floes"
    segF = morph_split_floes(segF_binarized, cloud_mask; p.floe_splitting_settings...)
    
    # Alternative: test whether the final segF_binarized is helpful
    # segF = morph_split_floes(ice_intersect, cloud_mask; max_fill_area=1, min_area_opening=20, opening_strel=se_disk(4))

    #### TBD: Test this subroutine and make into a function
    # Will need to remove small objects or merge them with prune_segments

    # Restore the original boundaries
    bw = .!(segF_binarized .|| segF) # Use the segF results to fill holes in segF binarized
    dist = .- distance_transform(feature_transform(bw))
    markers = label_components(segF)
    labels = labels_map(watershed(dist, markers))

    # bw2 = original k-means result ### TBD: Test which binarized image comes closest to the true floe boundaries
    bw2 = .!(kmeans_result .|| segF)
    labels[bw2] .= 0
    labels .= labels .* dilate(markers .> 0, se_disk(5))

    # add removal of too-small objects here. This needs to be label based not binary
    labels .= labels .* imfill(labels .> 0, (0, p.min_floe_size))

    # add removal of too-large objects here

    @info "Labeling floes"
    # labels = label_components(segF)
    labels = label_components(labels) # reset labels

    # Return the original truecolor image, segmented
    segments = SegmentedImage(truecolor, labels)

    if !isnothing(intermediate_results_callback)
        colorview_truecolor = view_seg(SegmentedImage(truecolor, labels))
        colorview_falsecolor = view_seg(SegmentedImage(falsecolor, labels))
        intermediate_results_callback(;
            truecolor,
            falsecolor,
            coastal_buffer_mask,
            cloud_mask,
            prelim_ice_mask,
            sharpened_grayscale_image=preproc_gray,
            bright_ice_mask=p.cluster_selection_algorithm(falsecolor_image),
            # ice_water_discrim=Gray.(ice_water_discrim), # TODO: Output of discriminate ice water can be forced to be Gray
            segA=kmeans_result,
            segAB_intersect=Gray.(ice_intersect),
            segB_enhanced_gray=Gray.(brightened_gray),
            watersheds_product=watersheds_product,
            segF_reconst_gray=morphed_grayscale,
            segF=segF_binarized, 
            final_floes = Gray.(labels_map(segments) .> 0),
            segment_mean_falsecolor=colorview_falsecolor,
            segment_mean_truecolor=colorview_truecolor,
            ) # Add figure that overlays the segments
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
function _reconstruct(sharpened_grayscale_image, landmask; strel=strel_diamond((5, 5)))
    markers = complement.(dilate(sharpened_grayscale_image, strel))
    mask = complement.(sharpened_grayscale_image)
    reconstructed_grayscale = mreconstruct(dilate, markers, mask)
    # dealing with complements! Ice is black in the constructed image, so we have to reverse it
    # to apply the landmask.
    reconstructed_grayscale = complement.(
        apply_landmask(
            complement.(reconstructed_grayscale),
             landmask)
    )
    return reconstructed_grayscale
end

"""
    clean_binary_floes(bw_img; min_opening_area=50)

Refine a binarized ice floe image (floes=white, leads/background/water=black) using morphological operations.

"""
function clean_binary_floes(bw_img; min_opening_area=50, strel=se_disk(3), min_object_size=50) # 
    # erode to separate objects
    out = erode(bw_img, strel)
    out .= imfill(out, (0, min_object_size)) 
    out .= closing(bw_img)
    out .= dilate(out, strel)
    out .= .!imfill(.!out, (0, min_opening_area))
    return out
end

"""
    clean_binary_floes_la2019

Version of post-segmentation cleanup used in LopezAcosta2019
"""
function clean_binary_floes_la2019(bw_img; min_opening_area=50)
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

""" # TODO: Sensitivity test. Does having the segB binarization in the workflow add value?
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
    img; # Only used to generate the segmented image for the output.
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
    stitched_labels = labels_map(stitch_clusters(wseg, tiles, min_overlap, grayscale_threshold))
    stitched_labels[.!binary_floes] .= 0
    return SegmentedImage(img, stitched_labels)
end

"""
    se_disk(r)

    Generate an approximately circular structuring element with radius r. For small r, this will be somewhat diamond-shaped.

""" # #TODO add simple example to docs, add to special strels, and in future, optimize the extreme filter for this shape
function se_disk(r)
    se = [sum(abs.(c.I .- (r + 1)) .^ 2) for c in CartesianIndices((2*r + 1, 2*r + 1))]
    return sqrt.(se) .<= r
end

function reconstruct_and_mask(
    grayscale_img::AbstractArray{<:Union{AbstractGray, TransparentGray}},
    watershed_boundary::BitMatrix,
    ice_intersect::BitMatrix,
    landmask::BitMatrix;
    se=strel_diamond((5,5)),
    min_area_opening=20
)
    ice_mask = .!watershed_boundary .* ice_intersect
    ice_mask .= .!area_opening(ice_mask; min_area=min_area_opening, connectivity=2)

    reconst_gray = dilate(grayscale_img, se)
    mreconstruct!(
        dilate, reconst_gray, complement.(reconst_gray), complement.(grayscale_img)
    )
    reconst_gray[ice_mask .== 0] .= 0 # better before or after reconstruction?

    # apply_landmask!(reconst_gray, landmask) # does this need to be here?
    return reconst_gray
end

"""
    morph_split_floes(binary_img; max_fill_area=1, min_area_opening=20, opening_strel=se_disk(4))

Separate floes in a binary image using hbreak, branch, imfill, bottom hat transform, and area opening.
Based on Lopez-Acosta et al. 2019, 2021.

"""
function morph_split_floes(binary_img, cloudmask; max_fill_area=1, min_area_opening=20, opening_strel=se_disk(4))
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

function clean_binary_floes2(binary_img, icemask, cloudmask; strel=strel_box((5,5)), max_fill=100)
    out = deepcopy(binary_img)
    eroded_img = erode(out, strel)
    filled = fill_holes(eroded_img, strel_diamond((3,3))) # Test how permissive this is.
    filled .= filled .&& (icemask .|| cloudmask)
    filled .= .!imfill(.!filled, (0, max_fill))
    out[filled .> 0] .= 1
    return out
end

end
