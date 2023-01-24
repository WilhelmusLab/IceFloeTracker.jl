"""
    segmentation_F(
    segmentation_C_ice_mask::BitMatrix,
    segmentation_B_not_ice_mask::BitMatrix,
    watershed_intersect::BitMatrix,
    cloudmask::BitMatrix,
    landmask::BitMatrix,
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
- `landmask.jl`: bitmatrix landmask for region of interest
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
    segmentation_B_not_ice_mask = IceFloeTracker.apply_landmask(
        segmentation_B_not_ice_mask, landmask
    )
    #blobs_applied
    ice_mask_watershed_applied = .!watershed_intersect .* segmentation_C_ice_mask
    #BW1
    ice_mask_watershed_opened = ImageMorphology.area_opening(
        ice_mask_watershed_applied; min_area=20, connectivity=2
    )
    #leads
    ice_leads = .!ice_mask_watershed_opened
    #Iobrd2
    not_ice_dilated = IceFloeTracker.MorphSE.dilate(
        segmentation_B_not_ice_mask; dims=IceFloeTracker.MorphSE.strel_diamond((3, 3))
    )
    #Iobrcbr2
    not_ice_reconstructed = IceFloeTracker.MorphSE.mreconstruct(
        IceFloeTracker.MorphSE.dilate,
        complement.(not_ice_dilated),
        complement.(segmentation_B_not_ice_mask)
    )
    #Z(he)
    reconstructed_leads = float64.(not_ice_reconstructed .* ice_leads) .+ (60 / 255)
    #BW_final
    leads_segmented = IceFloeTracker.kmeans_segmentation(reconstructed_leads, ice_labels)
    println("Done with k-means segmentation")
    #BW_final1_blobs_applied
    leads_segmented_watershed_applied = leads_segmented .* .!watershed_intersect
    #BW_final1  
    leads_segmented_broken =
        IceFloeTracker.hbreak(leads_segmented_watershed_applied) .* .!watershed_intersect
    #BW_final2
    leads_branched = IceFloeTracker.branch(leads_segmented_broken) .* .!watershed_intersect
    #BW_final3
    leads_filled = .!ImageMorphology.imfill(.!leads_branched, 0:1) .* .!watershed_intersect

    #BW_final4
    leads_opened = ImageMorphology.area_opening(leads_filled; min_area=20, connectivity=2)

    println("Done with area opening")
    #BW_final4(2)
    leads_opened_branched = IceFloeTracker.branch(leads_opened) 
    #BW_final4 
    leads_filled = IceFloeTracker.MorphSE.fill_holes(leads_opened_branched)
    #BW_final4_bothat
    leads_bothat = IceFloeTracker.MorphSE.bothat(
        leads_filled;
         dims=IceFloeTracker.MorphSE.strel_diamond((5,5))
    )
    #BW_final4(3)
    leads =
        convert(BitMatrix, (complement.(leads_bothat) .* leads_opened_branched)) .*
        .!watershed_intersect
    #BW1
    leads_bothat_opened = ImageMorphology.area_opening(leads; min_area=20, connectivity=2)
    #leads_bothat_opened = IceFloeTracker.branch(leads_bothat_opened)
    #BW2
    leads_bothat_filled = IceFloeTracker.MorphSE.fill_holes(leads_bothat_opened) .* cloudmask
    #BW2
    leads_bothat_masked = IceFloeTracker.branch(leads_bothat_filled ) #prune
    #BW3
    leads_cloudmasked_filled = IceFloeTracker.MorphSE.fill_holes(leads_bothat_masked)
    #BW4
    leads_masked_branched = IceFloeTracker.branch(leads_cloudmasked_filled)
    #BW5
    floes_erode = IceFloeTracker.MorphSE.erode(leads_masked_branched; dims=IceFloeTracker.se_disk4())
    floes_erode = IceFloeTracker.prune(IceFloeTracker.branch(floes_erode))
    #BW6
    floes_dilate = IceFloeTracker.MorphSE.dilate(floes_erode, IceFloeTracker.se_disk4())
    floes_opened = IceFloeTracker.prune(IceFloeTracker.branch(floes_dilate))
    floes_reconstructed = IceFloeTracker.MorphSE.mreconstruct(IceFloeTracker.MorphSE.dilate, leads_masked_branched, floes_dilate)

    return floes_dilate
end
