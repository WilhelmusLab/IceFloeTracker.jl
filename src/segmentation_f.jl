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
    segmentation_B_not_ice_mask::Matrix{Gray{Float64}},
    watershed_intersect::BitMatrix,
    cloudmask::BitMatrix,
    landmask::BitMatrix,
    ice_labels::Vector{Int64};
    lower_min_area_opening::Int64=20,
    upper_min_area_opening::Int64=200,
)::BitMatrix
    segmentation_B_not_ice_mask = IceFloeTracker.apply_landmask(segmentation_B_not_ice_mask, landmask)
#blobs_applied
    ice_mask_watershed_applied = .!watershed_intersect .* segmentation_C_ice_mask
#BW1
    ice_mask_watershed_opened = ImageMorphology.area_opening(
        ice_mask_watershed_applied, min_area=lower_min_area_opening
    )
#leads
    ice_leads = ifelse.(.!ice_mask_watershed_opened .== 0, 0.0, 1)
#Iobrd2
    not_ice_dilated = IceFloeTracker.MorphSE.dilate(
        segmentation_B_not_ice_mask, dims=IceFloeTracker.MorphSE.strel_diamond((5, 5))
    )
#Iobrcbr2
    not_ice_reconstructed = ImageMorphology.opening(
        complement.(not_ice_dilated), dims=complement.(segmentation_B_not_ice_mask)
    )
#Iobrcbr3
    reconstructed_leads = float64.(not_ice_reconstructed .* (ice_leads .+ (60 / 255)))
#BW_final
    leads_segmented, _, _ = IceFloeTracker.segmentation_A(
        reconstructed_leads, cloudmask, ice_labels
    )
    println("Done with k-means segmentation")
#BW_final1_blobs_applied
    leads_segmented_watershed_applied = leads_segmented .* .!watershed_intersect
#BW_final1  
    leads_segmented_broken =
    IceFloeTracker.hbreak!(leads_segmented_watershed_applied) .* .!watershed_intersect
#BW_final2
    leads_branched =
        IceFloeTracker.branch(leads_segmented_broken) .* .!watershed_intersect
#BW_final3
    leads_filled = ImageMorphology.imfill(.!leads_branched, 0:1) .* .!watershed_intersect
#BW_final4
    leads_opened =
        ImageMorphology.area_opening(.!leads_filled; min_area=lower_min_area_opening)
    println("Done with area opening")
#BW_final4(2)
    leads_opened_branched = IceFloeTracker.branch(leads_opened) #BW_final4 
#BW_final4_bothat
    leads_bothat = ImageMorphology.bothat(leads_opened_branched, dims=IceFloeTracker.MorphSE.strel_diamond((5, 5)))
#BW_final4(3)
    leads = convert(BitMatrix, (complement.(leads_bothat) .* leads_opened_branched)
    )
#BW1
    leads_bothat_opened =
        ImageMorphology.area_opening(leads, min_area=upper_min_area_opening)
#BW2
    leads_bothat_filled = 
        ImageMorphology.imfill(.!leads_bothat_opened, 0:upper_min_area_opening)
#BW2
    leads_bothat_masked = .!leads_bothat_filled .* cloudmask
#BW3
    leads_cloudmasked_filled = ImageMorphology.imfill(.!leads_bothat_masked, 0:upper_min_area_opening)
#BW4
    leads_masked_branched = IceFloeTracker.branch(.!leads_cloudmasked_filled)
#BW5
    floes_opened = ImageMorphology.opening(
        leads_masked_branched, dims=IceFloeTracker.se_disk4()
    )
    #isolated_floes = ImageMorphology.opening(leads_masked_branched, dims=floes_opened)
    return floes_opened
end
