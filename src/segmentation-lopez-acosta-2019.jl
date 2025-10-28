abstract type IceFloeSegmentationAlgorithm end

struct LopezAcosta2019 <: IceFloeSegmentationAlgorithm
    landmask_structuring_element::AbstractMatrix{Bool}
end

function LopezAcosta2019(; landmask_structuring_element=make_landmask_se())
    return LopezAcosta2019(landmask_structuring_element)
end

function (p::LopezAcosta2019)(
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
    cloudmask = create_cloudmask(falsecolor_image)

    # 2. Intermediate images
    @info "Finding ice labels"
    ice_mask = find_ice_mask(falsecolor_image, landmask_imgs.dilated)

    @info "Sharpening truecolor image"
    # a. apply imsharpen to truecolor image using non-dilated landmask
    sharpened_truecolor_image = imsharpen(truecolor_image, landmask_imgs.non_dilated)
    # b. apply imsharpen to sharpened truecolor img using dilated landmask
    sharpened_gray_truecolor_image = imsharpen_gray(
        sharpened_truecolor_image, landmask_imgs.dilated
    )

    @info "Normalizing truecolor image"
    normalized_image = normalize_image(
        sharpened_truecolor_image, sharpened_gray_truecolor_image, landmask_imgs.dilated
    )

    # Discriminate ice/water
    @info "Discriminating ice/water"
    ice_water_discrim = discriminate_ice_water(
        falsecolor_image, normalized_image, copy(landmask_imgs.dilated), cloudmask
    )

    # 3. Segmentation
    @info "Segmenting floes part 1/3"
    segA = segmentation_A(
        segmented_ice_cloudmasking(ice_water_discrim, cloudmask, ice_mask)
    )

    # segmentation_B
    @info "Segmenting floes part 2/3"
    segB = segmentation_B(sharpened_gray_truecolor_image, cloudmask, segA)

    # Process watershed in parallel using Folds
    @info "Building watersheds"
    # container_for_watersheds = [landmask_imgs.non_dilated, similar(landmask_imgs.non_dilated)]

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
            sharpened_gray_truecolor_image=sharpened_gray_truecolor_image,
            normalized_image=normalized_image,
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



Generates an image with ice floes apparent after filtering and combining previously processed versions of falsecolor and truecolor images from the same region of interest. Returns an image ready for segmentation to isolate floes.

Note: This function mutates the landmask object to avoid unnecessary memory allocation. If you need the original landmask, make a copy before passing it to this function. Example: `discriminate_ice_water(falsecolor_image, normalized_image, copy(landmask_bitmatrix), cloudmask_bitmatrix)`

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
    clouds_channel = IceFloeTracker.create_clouds_channel(
        cloudmask_bitmatrix, falsecolor_image
    )
    falsecolor_image_band7 = @view(channelview(falsecolor_image)[1, :, :])

    # first define all of the image variations
    image_clouds = IceFloeTracker.apply_landmask(clouds_channel, landmask_bitmatrix) # output during cloudmask apply, landmasked
    image_cloudless = IceFloeTracker.apply_landmask(
        falsecolor_image_band7, landmask_bitmatrix
    ) # channel 1 (band 7) from source falsecolor image, landmasked
    image_floes = IceFloeTracker.apply_landmask(falsecolor_image, landmask_bitmatrix) # source false color reflectance, landmasked
    image_floes_view = channelview(image_floes)

    floes_band_2 = @view(image_floes_view[2, :, :])
    floes_band_1 = @view(image_floes_view[3, :, :])

    # keep pixels greater than intensity 100 in bands 2 and 1
    floes_band_2_keep = floes_band_2[floes_band_2 .> floes_threshold]
    floes_band_1_keep = floes_band_1[floes_band_1 .> floes_threshold]

    _, floes_bin_counts = ImageContrastAdjustment.build_histogram(floes_band_2_keep, nbins)
    _, vals = Peaks.findmaxima(floes_bin_counts)

    differ = vals / (maximum(vals))
    proportional_intensity = sum(differ .> differ_threshold) / length(differ) # finds the proportional intensity of the peaks in the histogram

    # compute kurtosis, skewness, and standard deviation to use in threshold filtering
    kurt_band_2 = kurtosis(floes_band_2_keep)
    skew_band_2 = skewness(floes_band_2_keep)
    kurt_band_1 = kurtosis(floes_band_1_keep)
    standard_dev = stdmult(⋅, normalized_image)

    # find the ratio of clouds in the image to use in threshold filtering
    _, clouds_bin_counts = ImageContrastAdjustment.build_histogram(image_clouds .> 0)
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

    # reusing memory allocated in landmask_bitmatrix
    # used to be mask_image_clouds
    @. landmask_bitmatrix = (
        image_clouds < mask_clouds_lower || image_clouds > mask_clouds_upper
    )

    # reusing image_cloudless - used to be band7_masked
    @. image_cloudless = image_cloudless * !landmask_bitmatrix

    # reusing normalized_image_copy - used to be ice_water_discriminated_image
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

    segmented_ice_filled = IceFloeTracker.fill_holes(segmented_bridged)

    diff_matrix = segmented_ice_opened .!= segmented_ice_filled

    segmented_A = segmented_ice_cloudmasked .|| diff_matrix

    return segmented_A
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
    gray_image::Matrix{Gray{Float64}},
    cloudmask::BitMatrix,
    ice_labels::Union{Vector{Int64},BitMatrix,AbstractArray{<:Gray}},
)::BitMatrix
    segmented_ice = IceFloeTracker.kmeans_segmentation(gray_image, ice_labels)
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
    features = Images.feature_transform(.!intermediate_segmentation_image)
    distances = 1 .- Images.distance_transform(features)
    seg_mask = ImageSegmentation.hmin_transform(distances, 2)
    seg_mask_bool = seg_mask .> 0
    markers = Images.label_components(seg_mask_bool)
    segment = ImageSegmentation.watershed(distances, markers)
    labels = ImageSegmentation.labels_map(segment)
    borders = Images.isboundary(labels)
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
