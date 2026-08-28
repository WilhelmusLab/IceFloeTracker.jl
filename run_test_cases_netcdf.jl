using IceFloeTracker
using NCDatasets
using Images

##### Notes #####
# The test cases saved the results into netcdf files, so I can check
#  whether the issues from the FSPipeline version on Main is corrected with the newer filtering functions.

ds = NCDataset("/Users/dmw/Downloads/case_studies/beaufort_sea-large.250m.2008-06-08.aqua.nc");

# Extract the true color and false color images
tc_img_data = ds["modis_truecolor"][:,:,:] ./ 255;
fc_img_data = ds["modis_falsecolor"][:,:,:] ./ 255;

# Need to reformat it to an image
tc_img = colorview(RGBA, [transpose(tc_img_data[:, :, i]) for i in 1:4]...);
fc_img = colorview(RGBA, [transpose(fc_img_data[:, :, i]) for i in 1:4]...);
lm_data = transpose(ds["landmask"][:, :]);
landmask = coalesce.(lm_data, 0) .> 0

cm_data = transpose(ds["coastal_buffer_mask"][:, :]);
coast_buffer = coalesce.(cm_data, 0) .> 0


# Now that the data is loaded and formatted, we can test the components of the FSPipeline
# This is a larger image than the Greenland Sea cases -- 4224 by 6000 pixels. It also contains larger
# floes than that case.
# Default settings, 400 pixel tiles: 5,597 detected objects, 238 seconds on my Mac Mini
# Changing the size to 1000 and dropping the probability-based filter: 279 seconds, 6,179 shapes. However, most of these are tiny and likely are artifacts.
floe_filtering_params = (
    min_floe_size=100,
    max_floe_size=90_000,
    boundary_radius=15,
    min_circularity=0.3,
    min_solidity=0.7,
    min_reflectance=0.4,
    min_contrast=0.01, # should this be "boundary contrast"?
    min_probability=0.,
)
@time begin
    segment = IceFloeTracker.FSPipeline.Segment(;tile_size_pixels=1000, floe_filtering_params=floe_filtering_params)
    segmented_image = segment(tc_img, fc_img, landmask);
end

save("/Users/dmw/Downloads/case_studies/beaufort_sea-large.250m.2008-06-08.aqua.fs_pipeline_update_tilesize_1000px.png", view_seg_random(segmented_image))

# Result with previous version: 3,887 floes.
labeled_data = Int64.(coalesce.(transpose(ds["labeled_image"][:, :]), 0));
segmented_initial = SegmentedImage(tc_img, labeled_data)
save("/Users/dmw/Downloads/case_studies/beaufort_sea-large.250m.2008-06-08.aqua.fs_pipeline_current.png",
 view_seg_random(segmented_initial))

# The LopezAcosta2019Tiling method at first glance is outperforming the FSPipeline for retrieved floes, while
# taking longer to run. 602 seconds and 15,484 floes, though most of these floes are likely false positives.
@time begin
    segment = IceFloeTracker.LopezAcosta2019Tiling.Segment()
    segmented_image = segment(tc_img, fc_img, landmask);
end

save("/Users/dmw/Downloads/case_studies/beaufort_sea-large.250m.2008-06-08.aqua.LopezAcosta2019Tiling.png", view_seg_random(segmented_image))

#### Some questions:
# How consistently is this method better than mine? Is the issue the choice of tile size, or the filtering step?
# Is the issue the calibration against small images? Perhaps my validation method is all wrong.
labeled_split = FSPipeline.dist_morph_split(labels_map(segmented_image) .> 0; max_hole_fill=2000, max_distance=5, max_expand=3);
save("/Users/dmw/Downloads/case_studies/beaufort_sea-large.250m.2008-06-08.aqua.LopezAcosta2019Tiling_FSSplit.png", 
    view_seg_random(SegmentedImage(tc_img, labeled_split)))

filtered_floes = FSPipeline.filter_floes(labeled_split, coast_buffer, Watkins2026CloudMask()(fc_img), fc_img)
FSPipeline.keep_labels!(labeled_split, filtered_floes.label)   
save("/Users/dmw/Downloads/case_studies/beaufort_sea-large.250m.2008-06-08.aqua.LopezAcosta2019Tiling_FSSplitFiltered.png", 
    view_seg_random(SegmentedImage(tc_img, labeled_split)))



#### Neither of these looks as good as I want it to look.
# Let's check what stage things go off the rails.

cloud_mask = Watkins2026CloudMask()(fc_img)
truecolor_image = float64.(RGB.(tc_img))
falsecolor_image = float64.(RGB.(fc_img))
land_mask = landmask .> 0 # make sure it's a bitmatrix, in case it's passed as Gray
apply_landmask!(tc_img, land_mask)
apply_landmask!(fc_img, land_mask)

n, m = size(truecolor_image)
tile_size_pixels = 1200
nmin, nmax = extrema([n, m])
tile_size_pixels > nmax && begin
    @warn "Tile size too large; clamping to min(height, width)."
    tile_size_pixels = nmin
end

(nr, nc) = round.(Int, size(tc_img) ./ tile_size_pixels)
tiles = get_tiles(tc_img; rblocks=nr, cblocks=nc)

# 2. Intermediate images - apply coastal buffer and cloud mask
joint_mask = coast_buffer .|| cloud_mask
tc_masked = apply_landmask(tc_img, joint_mask)
fc_masked = apply_landmask(fc_img, joint_mask)

# First check for sufficient non-land and non-cloud pixels
filtered_tiles = filter(
    t -> sum(.!joint_mask[t...]) > 300, tiles
);

# Then check for sufficient possible sea ice pixels
prelim_ice_mask = IceDetectionBrightnessMidpoint(; minimum_reflectance=0.3)(Gray.(red.(tc_masked)), filtered_tiles);
filtered_tiles = filter(
    t -> sum(prelim_ice_mask[t...]) > 300, filtered_tiles
);
water_mask = .!(prelim_ice_mask .|| cloud_mask .|| land_mask);

preproc_gray = float64.(
    FSPipeline.preprocessing_algorithm(truecolor_image, land_mask, filtered_tiles)
);

mosaicview(tc_img, Gray.(prelim_ice_mask), Gray.(water_mask), nrow=1)

adaptive_result = apply_landmask(
        binarize(preproc_gray, AdaptiveThreshold(; window_size=800, percentage=0)) .> 0,
        water_mask .|| land_mask
);

kmeans_result = kmeans_binarization(
    preproc_gray, fc_masked, filtered_tiles; FSPipeline.kmeans_params...
);

mosaicview(tc_img, #view_seg_random(segmented_initial),
 Gray.(adaptive_result), Gray.(kmeans_result), ncol=2)
clean_and_split(r) = FSPipeline.dist_morph_split(
        FSPipeline.clean_binary_floes(
            r, prelim_ice_mask, cloud_mask; FSPipeline.cleanup_binary_params...
        ); max_hole_fill = 2000, max_distance = 15, max_expand = 3
)
labeled_images = clean_and_split.([kmeans_result, adaptive_result]);
mosaicview(view_seg_random(SegmentedImage(tc_img, labeled_images[1])),
    view_seg_random(SegmentedImage(tc_img, labeled_images[2])), ncol=1)

# tiled watershed