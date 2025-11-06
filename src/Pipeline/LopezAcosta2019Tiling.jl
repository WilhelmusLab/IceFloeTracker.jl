
module LopezAcosta2019Tiling

import Images:
    Images,
    area_opening,
    watershed,
    imfilter,
    isboundary,
    distance_transform,
    feature_transform,
    label_components,
    strel_diamond,
    strel_box,
    AbstractRGB,
    TransparentRGB,
    AbstractGray,
    dilate,
    mreconstruct,
    RGB,
    Gray,
    float64,
    red,
    green,
    blue,
    adjust_histogram,
    GammaCorrection,
    opening,
    centered,
    labels_map,
    local_maxima,
    SegmentedImage,
    segment_mean
    

import ..skimage: sk_morphology
import ..ImageUtils: get_brighten_mask, to_uint8, imcomplement, imbrighten, get_tiles
import ..Filtering: histeq, unsharp_mask, conditional_histeq, rgb2gray
import ..Preprocessing:
    apply_landmask,
    apply_landmask!,
    apply_cloudmask,
    create_cloudmask,
    create_landmask,
    LopezAcostaCloudMask
import ..Morphology:
    hbreak!,
    hbreak,
    fill_holes,
    morph_fill,
    reconstruct,
    branch,
    bridge,
    se_disk4,
    se_disk2,
    se_disk20,
    imextendedmin,
    impose_minima,
    imregionalmin
import ..Segmentation:
    IceFloeSegmentationAlgorithm,
    tiled_adaptive_binarization,
    kmeans_segmentation,
    get_ice_masks

# Sample input parameters expected by the main function
cloud_mask_thresholds = (
    prelim_threshold=110.0 / 255.0,
    band_7_threshold=200.0 / 255.0,
    band_2_threshold=190.0 / 255.0,
    ratio_lower=0.0,
    ratio_offset=0.0,
    ratio_upper=0.75,
)

adapthisteq_params = (
    white_threshold=25.5, entropy_threshold=4, white_fraction_threshold=0.4
)

adjust_gamma_params = (gamma=1.5, gamma_factor=1.3, gamma_threshold=220)

structuring_elements = (
    se_disk1=collect(strel_diamond((3, 3))), se_disk2=se_disk2(), se_disk4=se_disk4()
)

unsharp_mask_params = (radius=10, amount=2.0, factor=255.0)

brighten_factor = 0.1

ice_masks_params = (
    band_7_threshold=5 / 255,
    band_2_threshold=230 / 255,
    band_1_threshold=240 / 255,
    band_7_threshold_relaxed=10 / 255,
    band_1_threshold_relaxed=190 / 255,
    possible_ice_threshold=75 / 255,
    k=3, # number of clusters for kmeans segmentation
)

prelim_icemask_params = (radius=10, amount=2, factor=0.5)

@kwdef struct Segment <: IceFloeSegmentationAlgorithm
    tile_settings = (; rblocks=2, cblocks=2)
    cloud_mask_thresholds = cloud_mask_thresholds
    adapthisteq_params = adapthisteq_params
    adjust_gamma_params = adjust_gamma_params
    structuring_elements = structuring_elements
    unsharp_mask_params = unsharp_mask_params
    ice_masks_params = ice_masks_params
    prelim_icemask_params = prelim_icemask_params
    brighten_factor = brighten_factor
end

function (p::Segment)(
    truecolor::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    falsecolor::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    landmask::AbstractArray{<:Union{AbstractGray,AbstractRGB,TransparentRGB}};
    intermediate_results_callback::Union{Nothing,Function}=nothing,
)
    _landmask = create_landmask(landmask, strel_box((3, 3))) # smaller strel than in some test cases
    tiles = get_tiles(truecolor; p.tile_settings...)

    ref_image = RGB.(falsecolor)  # TODO: remove this typecast
    true_color_image = RGB.(truecolor)  # TODO: remove this typecast

    begin
        @debug "Step 1/2: Create and apply cloudmask to reference image"

        cloudmask = create_cloudmask(
            ref_image, LopezAcostaCloudMask(cloud_mask_thresholds...)
        )
        ref_img_cloudmasked = apply_cloudmask(ref_image, cloudmask)
    end

    begin
        @debug "Step 3: Tiled adaptive histogram equalization"
        clouds_red = to_uint8(float64.(red.(ref_img_cloudmasked) .* 255))
        clouds_red[_landmask.dilated] .= 0

        rgbchannels = conditional_histeq(
            true_color_image, clouds_red, tiles; adapthisteq_params...
        )

        gammagreen = @view rgbchannels[:, :, 2]
        equalized_gray = rgb2gray(rgbchannels)
    end

    begin
        @debug "Step 4: Remove clouds from equalized_gray"
        masks = [f.(ref_img_cloudmasked) .== 0 for f in [red, green, blue]]
        combo_mask = reduce((a, b) -> a .& b, masks)
        equalized_gray[combo_mask] .= 0
    end

    begin
        @debug "Step 5: unsharp_mask on equalized_gray and reconstruct"
        sharpened = to_uint8(unsharp_mask(equalized_gray, unsharp_mask_params...))
        equalized_gray_sharpened_reconstructed = reconstruct(
            sharpened, structuring_elements.se_disk1, "dilation", true
        )
        equalized_gray_sharpened_reconstructed[_landmask.dilated] .= 0
    end

    # TODO: Steps 6 and 7 can be done in parallel as they are independent
    begin
        @debug "Step 6: Repeat step 5 with equalized_gray (landmasking, no sharpening)"
        equalized_gray_reconstructed = deepcopy(equalized_gray)
        apply_landmask!(equalized_gray_reconstructed, _landmask.dilated)

        equalized_gray_reconstructed = reconstruct(
            equalized_gray_reconstructed, structuring_elements.se_disk1, "dilation", true
        )
        apply_landmask!(equalized_gray_reconstructed, _landmask.dilated)
    end

    begin
        @debug "Step 7: Brighten equalized_gray"
        brighten = get_brighten_mask(equalized_gray_reconstructed, gammagreen)
        apply_landmask!(equalized_gray, _landmask.dilated)
        equalized_gray .= imbrighten(equalized_gray, brighten, brighten_factor)
    end

    begin
        @debug "Step 8: Get morphed_residue and adjust its gamma"
        morphed_residue = clamp.(equalized_gray - equalized_gray_reconstructed, 0, 255)

        agp = adjust_gamma_params
        equalized_gray_sharpened_reconstructed_adjusted = imcomplement(
            adjustgamma(equalized_gray_sharpened_reconstructed, agp.gamma)
        )
        adjusting_mask =
            equalized_gray_sharpened_reconstructed_adjusted .> agp.gamma_threshold
        morphed_residue[adjusting_mask] .=
            to_uint8.(morphed_residue[adjusting_mask] .* agp.gamma_factor)
    end

    begin
        @debug "Step 9: Get preliminary ice masks"
        binarized_tiling =
            tiled_adaptive_binarization(
                Gray.(morphed_residue ./ 255),
                tiles;
                minimum_window_size=32,
                threshold_percentage=15,
            ) .> 0


        
        prelim_icemask = get_ice_masks(
            ref_image,
            Gray.(morphed_residue / 255),
            _landmask.dilated,
            tiles;
            ice_masks_params...,
        )
    end

    begin
        @debug "Step 10: Get segmentation mask from preliminary icemask"
        # Fill holes function in get_segment_mask a bit more aggressive than Matlabs
        segment_mask = get_segment_mask(prelim_icemask, binarized_tiling)
    end

    begin
        @debug "Step 11: Get local_maxima_mask and L0mask via watershed"
        local_maxima_mask, L0mask = watershed2(
            morphed_residue, segment_mask, prelim_icemask
        )
    end

    begin
        @debug "Step 12: Build icemask from all others"
        local_maxima_mask = to_uint8(local_maxima_mask * 255)
        prelim_icemask2 = _regularize(
            morphed_residue,
            local_maxima_mask,
            segment_mask,
            L0mask,
            structuring_elements.se_disk1;
            prelim_icemask_params...,
        )
    end

    begin
        @debug "Step 13: Get improved icemask"
        icemask = get_ice_masks(
            ref_image,
            Gray.(prelim_icemask2 ./ 255),
            _landmask.dilated,
            tiles;
            ice_masks_params...,
        )
    end

    begin
        @debug "Step 14: Get final mask"
        se = structuring_elements
        se_erosion = se.se_disk1
        se_dilation = se.se_disk2
        binary_floe_masks = get_final(icemask, segment_mask, se_erosion, se_dilation)
    end

    begin
        @debug "Step 15: Converting mask to segmented image"
        labels = label_components(binary_floe_masks)
        segmented = SegmentedImage(truecolor, labels)
    end

    if !isnothing(intermediate_results_callback)
        segments_truecolor = SegmentedImage(truecolor, labels)
        segments_falsecolor = SegmentedImage(falsecolor, labels)
        intermediate_results_callback(;
            falsecolor,
            truecolor,
            ref_img_cloudmasked,
            gammagreen,
            equalized_gray,
            equalized_gray_sharpened_reconstructed,
            equalized_gray_reconstructed,
            morphed_residue,
            prelim_icemask,
            binarized_tiling,
            segment_mask,
            local_maxima_mask,
            L0mask,
            prelim_icemask2,
            icemask,
            binary_floe_masks,
            labels,
            segmented,
            segment_mean_truecolor=map(i -> segment_mean(segments_truecolor, i), labels),
            segment_mean_falsecolor=map(i -> segment_mean(segments_falsecolor, i), labels),
        )
    end

    return segmented
end




function get_holes(img, min_opening_area=20, se=se_disk4())
    _img = area_opening(img; min_area=min_opening_area)
    hbreak!(_img)

    out = branchbridge(_img)
    out = opening(out, centered(se))
    out = fill_holes(out)

    return out .!= _img
end

function fillholes!(img)
    img[get_holes(img)] .= true
    return nothing
end

function get_segment_mask(ice_mask, tiled_binmask)
    # TODO: Threads.@threads # sometimes crashes (too much memory?)
    for img in (ice_mask, tiled_binmask)
        fillholes!(img)
        img .= watershed1(img)
    end
    segment_mask = ice_mask .&& tiled_binmask
    return segment_mask
end

function branchbridge(img)
    img = branch(img)
    img = bridge(img)
    return img
end

function watershed1(bw::T) where {T<:Union{BitMatrix,AbstractMatrix{Bool}}}
    seg = -bwdist(.!bw)
    mask2 = imextendedmin(seg)
    seg = impose_minima(seg, mask2)
    cc = label_components(imregionalmin(seg), trues(3, 3))
    w = watershed(seg, cc)
    lmap = labels_map(w)
    return isboundary(lmap) .> 0
    #dmw: isboundary returns a thick boundary, whereas matlab uses a 1-pixel boundary.
end

"""
    bwdist(bwimg)

Distance transform for binary image `bwdist`.
"""
function bwdist(bwimg::AbstractArray{Bool})::AbstractArray{Float64}
    return distance_transform(feature_transform(bwimg))
end

function _reconst_watershed(morph_residue::Matrix{<:Integer}, se::Matrix{Bool}=se_disk20())
    mr_reconst = to_uint8(reconstruct(morph_residue, se, "erosion", false))
    mr_reconst .= to_uint8(reconstruct(mr_reconst, se, "dilation", true))
    mr_reconst .= imcomplement(mr_reconst)
    return mr_reconst
end

function watershed2(morph_residue, segment_mask, ice_mask)
    # TODO: reconfigure to use async tasks or threads
    # Task 1: Reconstruct morph_residue
    # task1 = Threads.@spawn begin
    mr_reconst = _reconst_watershed(morph_residue)
    mr_reconst = local_maxima(mr_reconst; connectivity=2) .> 0
    # end

    # Task 2: Calculate gradient magnitude
    # task2 = Threads.@spawn begin
    gmag = imgradientmag(histeq(morph_residue))
    # end

    # Wait for both tasks to complete
    # mr_reconst = fetch(task1)
    # gmag = fetch(task2)

    minimamarkers = mr_reconst .| segment_mask .| ice_mask
    gmag .= impose_minima(gmag, minimamarkers)
    cc = label_components(imregionalmin(gmag), trues(3, 3))
    w = watershed(morph_residue, cc)
    lmap = labels_map(w)
    return (fgm=mr_reconst, L0mask=isboundary(lmap) .> 0)
end

"""
    imgradientmag(img)

Compute the gradient magnitude of an image using the Sobel operator.
"""
function imgradientmag(img)
    h = centered([-1 0 1; -2 0 2; -1 0 1]')
    Gx_future = Threads.@spawn imfilter(img, h', "replicate")
    Gy_future = Threads.@spawn imfilter(img, h, "replicate")
    Gx = fetch(Gx_future)
    Gy = fetch(Gy_future)
    return hypot.(Gx, Gy)
end

"""
    regularize_fill_holes(img, local_maxima_mask, factor, segment_mask, L0mask)

Regularize `img` by:
    1. increasing the maxima of `img` by a factor of `factor`
    2. filtering `img` at positions where either `segment_mask` or `L0mask` are true
    3. filling holes

# Arguments
- `img`: The morphological residue image.
- `local_maxima_mask`: The local maxima mask.
- `factor`: The factor to apply to the local maxima mask.
- `segment_mask`: The segment mask -- intersection of bw1 and bw2 in first tiled workflow of `master.m`.
- `L0mask`: zero-labeled pixels from watershed.
"""
function regularize_fill_holes(img, local_maxima_mask, segment_mask, L0mask, factor)
    new2 = to_uint8(img .+ local_maxima_mask .* factor)
    new2[segment_mask .|| L0mask] .= 0
    return fill_holes(new2)
end

"""
    regularize_sharpening(img, L0mask, radius, amount, local_maxima_mask, factor, segment_mask, se)

Regularize `img` via sharpening, filtering, reconstruction, and maxima elevating.

# Arguments
- `img`: The input image.
- `L0mask`: zero-labeled pixels from watershed.
- `radius`: The radius of the unsharp mask.
- `amount`: The amount of unsharp mask.
- `local_maxima_mask`: The local maxima mask.
- `factor`: The factor to apply to the local maxima mask.
- `segment_mask`: The segment mask -- intersection of bw1 and bw2 in first tiled workflow of `master.m`.
"""
function regularize_sharpening(
    img, L0mask, local_maxima_mask, segment_mask, se, radius, amount, factor
)
    new3 = unsharp_mask(img, radius, amount, 255)
    new3[L0mask] .= 0
    new3 = reconstruct(new3, se, "dilation", false)
    new3[segment_mask] .= 0
    return to_uint8(new3 + local_maxima_mask .* factor)
end

function _regularize(
    morph_residue, local_maxima_mask, segment_mask, L0mask, se; factor, radius, amount
)
    reg_fill_holes = regularize_fill_holes(
        morph_residue, local_maxima_mask, segment_mask, L0mask, factor[1]
    )
    reg_sharpened = regularize_sharpening(
        reg_fill_holes,
        L0mask,
        local_maxima_mask,
        segment_mask,
        se,
        radius,
        amount,
        factor[end],
    )
    return reg_sharpened
end

"""
    get_final(img, label, segment_mask, se_erosion, se_dilation)

Final processing following the tiling workflow.

# Arguments
- `img`: The input image.
- `label`: Mode of most common label in the find_ice_labels workflow.
- `segment_mask`: The segment mask.
- `se_erosion`: structuring element for erosion.
- `se_dilation`: structuring element for dilation.
- `apply_segment_mask=true`: Whether to filter `img` the segment mask.

"""
function get_final(
    # this function is used in preprocessing_tiling
    img,
    segment_mask,
    se_erosion,
    se_dilation,
    apply_segment_mask::Bool=true,
)
    _img = hbreak(img)

    # slow for big images
    _img .= morph_fill(_img)

    # TODO: decide on criteria for applying segment mask
    apply_segment_mask && (_img[segment_mask] .= false)

    # tends to fill more than matlabs imfill
    _img .= fill_holes(_img)

    # marker image
    _img .= branch(_img)

    #= opening to remove noise while preserving shape/size
    Note the different structuring elements for erosion and dilation =#
    mask = sk_morphology.erosion(_img, se_erosion)
    mask .= sk_morphology.dilation(mask, se_dilation)

    # Restore shape of floes based on the cleaned up `mask`
    final = mreconstruct(dilate, _img, mask)
    return BitMatrix(final)
end

# TODO: Remove once the workflow is all normed images
function adjustgamma(img, gamma=1.5, asuint8=true)
    if maximum(img) > 1
        img = img ./ 255
    end

    adjusted = adjust_histogram(img, GammaCorrection(gamma))

    if asuint8
        adjusted = Int.(round.(adjusted * 255, RoundNearestTiesAway))
    end

    return adjusted
end

end
