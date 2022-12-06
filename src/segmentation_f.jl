"""
    segmentation_F(
    segmentation_C_ice_mask::BitMatrix,
    segmentation_B_not_ice_mask::BitMatrix,
    watershed_intersect::BitMatrix,
    cloudmask::BitMatrix,
    ice_labels::Vector{Int64};
    lower_min_area_opening::Int64=20,
    upper_min_area_opening::Int64=150,
)::BitMatrix

Cleans up past segmentation images with morphological operations, and applies the results of prior watershed segmentation, returning the final cleaned image for tracking with ice floes segmented and isolated. 

# Arguments
- `segmentation_C_ice_mask`: binary cloudmasked and landmasked intermediate file from segmentation_C (`segmented_ice`)
- `segmentation_B_not_ice_mask`: binary mask output from `segmentation_b.jl`
- `watershed_intersect`: ice pixels, output from `segmentation_d_e.jl` 
- `cloudmask.jl`: bitmatrix cloudmask for region of interest
- `ice_labels`: vector of pixel coordinates output from `find_ice_labels.jl`
- `lower_min_area_opening`: threshold used for area opening; pixel groups greater than threshold are retained
- `upper_min_area_opening`: threshold used for area opening; pixel groups greater than threshold are retained

"""
function segmentation_F(
    segmentation_C_ice_mask::BitMatrix,
    segmentation_B_not_ice_mask::BitMatrix,
    watershed_intersect::BitMatrix,
    cloudmask::BitMatrix,
    ice_labels::Vector{Int64};
    lower_min_area_opening::Int64=20,
    upper_min_area_opening::Int64=150,
)::BitMatrix
    ice_mask_watershed_applied = .!watershed_intersect .* segmentation_C_ice_mask
    ice_mask_watershed_opened = ImageMorphology.area_opening(
        ice_mask_watershed_applied; min_area=lower_min_area_opening
    )
    ice_leads = ifelse.(.!ice_mask_watershed_opened .== 0, 0.0, 1 + (60 / 255))
    not_ice_dilated = IceFloeTracker.MorphSE.dilate(
        segmentation_B_not_ice_mask; dims=IceFloeTracker.MorphSE.strel_diamond((3, 3))
    )
    not_ice_reconstructed = ImageMorphology.opening(
        complement.(not_ice_dilated); dims=complement.(segmentation_B_not_ice_mask)
    )
    reconstructed_leads = float64.(not_ice_reconstructed .* ice_leads)
    leads_segmented, _, _ = IceFloeTracker.segmentation_A(
        Gray.(reconstructed_leads), cloudmask, ice_labels
    )
    println("Done with k-means segmentation")
    leads_segmented_watershed_applied = leads_segmented .* .!watershed_intersect
    IceFloeTracker.hbreak!(leads_segmented_watershed_applied)
    leads_segmented_watershed_applied =
        leads_segmented_watershed_applied .* .!watershed_intersect
    println(typeof(leads_segmented_watershed_applied))
    leads_branched =
        IceFloeTracker.branch(leads_segmented_watershed_applied) .* .!watershed_intersect
    leads_filled = ImageMorphology.imfill(.!leads_branched, 0:5) .* .!watershed_intersect
    leads_opened = IceFloeTracker.branch(
        ImageMorphology.area_opening(.!leads_filled; min_area=upper_min_area_opening)
    )
    println("Done with area opening")
    leads_bothat = convert(
        BitMatrix, (complement.(ImageMorphology.bothat(.!leads_opened)) .* leads_opened)
    )
    leads_bothat_opened =
        ImageMorphology.area_opening(.!leads_bothat; min_area=lower_min_area_opening) .*
        cloudmask
    leads_bothat_filled = IceFloeTracker.branch(
        ImageMorphology.imfill(.!leads_bothat_opened, 0:10)
    )
    floes_opened = ImageMorphology.opening(
        leads_bothat_filled; dims=IceFloeTracker.se_disk4()
    )
    isolated_floes = ImageMorphology.opening(leads_bothat_filled; dims=floes_opened)
    return isolated_floes
end
