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
    ice_labels::Union{Vector{Int64},BitMatrix},
    cloudmask::BitMatrix,
    landmask::BitMatrix;
    min_area_opening::Int64=20,
)::BitMatrix
    IceFloeTracker.apply_landmask!(segmentation_B_not_ice_mask, landmask)

    ice_leads = .!segmentation_B_watershed_intersect .* segmentation_B_ice_intersect

    ice_leads .=
        .!ImageMorphology.area_opening(ice_leads; min_area=min_area_opening, connectivity=2)

    not_ice = IceFloeTracker.dilate(
        segmentation_B_not_ice_mask, IceFloeTracker.strel_diamond((5, 5))
    )

    IceFloeTracker.mreconstruct!(
        IceFloeTracker.dilate,
        not_ice,
        complement.(not_ice),
        complement.(segmentation_B_not_ice_mask),
    )

    reconstructed_leads = (not_ice .* ice_leads) .+ (60 / 255)

    leads_segmented =
        IceFloeTracker.kmeans_segmentation(reconstructed_leads, ice_labels) .*
        .!segmentation_B_watershed_intersect
    @info("Done with k-means segmentation")
    leads_segmented_broken = IceFloeTracker.hbreak(leads_segmented)

    leads_branched = IceFloeTracker.branch(leads_segmented_broken)

    leads_filled = .!ImageMorphology.imfill(.!leads_branched, 0:1)

    leads_opened = IceFloeTracker.branch(
        ImageMorphology.area_opening(
            leads_filled; min_area=min_area_opening, connectivity=2
        ),
    )

    leads_bothat =
        IceFloeTracker.bothat(leads_opened, IceFloeTracker.strel_diamond((5, 5))) .> 0.499

    leads = convert(BitMatrix, (complement.(leads_bothat) .* leads_opened))

    ImageMorphology.area_opening!(leads, leads; min_area=min_area_opening, connectivity=2)

    # dmw: replace multiplication with apply_cloudmask
    leads_bothat_filled = (IceFloeTracker.fill_holes(leads) .* .!cloudmask)
    # leads_bothat_filled = apply_cloudmask(IceFloeTracker.fill_holes(leads), cloudmask)
    floes = IceFloeTracker.branch(leads_bothat_filled)

    floes_opened = IceFloeTracker.opening(floes, centered(IceFloeTracker.se_disk4()))

    IceFloeTracker.mreconstruct!(IceFloeTracker.dilate, floes_opened, floes, floes_opened)

    return floes_opened
end
