"""
    segmentation_F(
    segmentation_C_ice_mask::BitMatrix,
    segmentation_B_not_ice_mask::BitMatrix,
    watershed_intersect::BitMatrix,
    cloudmask::BitMatrix,
    landmask::BitMatrix,
    ice_labels::Vector{Int64};
    min_area_opening::Int64=20,
)::BitMatrix

Cleans up past segmentation images with morphological operations, and applies the results of prior watershed segmentation, returning the final cleaned image for tracking with ice floes segmented and isolated. 

# Arguments
- `segmentation_C_ice_mask`: binary cloudmasked, landmasked intermediate file from segmentation_C (`segmented_ice`)
- `segmentation_B_not_ice_mask`: binary mask output from `segmentation_b.jl`
- `watershed_intersect`: ice pixels, output from `segmentation_d_e.jl` 
- `cloudmask.jl`: bitmatrix cloudmask for region of interest
- `landmask.jl`: bitmatrix landmask for region of interest
- `ice_labels`: vector of pixel coordinates output from `find_ice_labels.jl`
- `min_area_opening`: threshold used for area opening; pixel groups greater than threshold are retained

"""
function segmentation_F(
    segmentation_C_ice_mask::BitMatrix,
    segmentation_B_not_ice_mask::Matrix{Gray{Float64}},
    watershed_intersect::BitMatrix,
    cloudmask::BitMatrix,
    landmask::BitMatrix,
    ice_labels::Vector{Int64};
    min_area_opening::Int64=20,
)::BitMatrix
    segmentation_B_not_ice_mask = IceFloeTracker.apply_landmask!(
        segmentation_B_not_ice_mask, landmask
    )

    ice_mask_watershed = .!watershed_intersect .* segmentation_C_ice_mask

    ice_mask_watershed = ImageMorphology.area_opening(
        ice_mask_watershed; min_area=min_area_opening, connectivity=2
    )

    ice_leads = .!ice_mask_watershed

    not_ice_dilated = IceFloeTracker.MorphSE.dilate(
        segmentation_B_not_ice_mask, IceFloeTracker.MorphSE.strel_diamond((5, 5))
    )

    not_ice_reconstructed = IceFloeTracker.MorphSE.mreconstruct(
        IceFloeTracker.MorphSE.dilate,
        complement.(not_ice_dilated),
        complement.(segmentation_B_not_ice_mask),
    )

    reconstructed_leads = (float64.(not_ice_reconstructed .* ice_leads) .+ (60 / 255))

    leads_segmented =
        IceFloeTracker.kmeans_segmentation(reconstructed_leads, ice_labels) .*
        .!watershed_intersect
    println("Done with k-means segmentation")
    leads_segmented_broken = IceFloeTracker.hbreak(leads_segmented)

    leads_branched = IceFloeTracker.branch(leads_segmented_broken)

    leads_filled = .!ImageMorphology.imfill(.!leads_branched, 0:1)

    leads_opened = ImageMorphology.area_opening(
        leads_filled; min_area=min_area_opening, connectivity=2
    )

    leads_opened_branched = IceFloeTracker.branch(leads_opened)

    leads_bothat =
        IceFloeTracker.MorphSE.bothat(
            leads_opened_branched, IceFloeTracker.MorphSE.strel_diamond((5, 5))
        ) .> 0.499

    leads = convert(BitMatrix, (complement.(leads_bothat) .* leads_opened_branched))

    leads_bothat_opened = ImageMorphology.area_opening(
        leads; min_area=min_area_opening, connectivity=2
    )

    leads_bothat_filled = (
        IceFloeTracker.MorphSE.fill_holes(leads_bothat_opened) .* cloudmask
    )

    leads_masked_branched = IceFloeTracker.branch(leads_bothat_filled)

    floes_opened = IceFloeTracker.MorphSE.opening(
        leads_masked_branched, IceFloeTracker.se_disk4()
    )

    floes_reconstructed = IceFloeTracker.MorphSE.mreconstruct(
        IceFloeTracker.MorphSE.dilate, leads_masked_branched, floes_opened
    )

    return floes_reconstructed
end
