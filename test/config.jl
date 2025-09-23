# Setting things up

## locate some files for the tests
test_data_dir = "./test_inputs"
test_output_dir = "./test_outputs"
truecolor_test_image_file = "$(test_data_dir)/beaufort-chukchi-seas_truecolor.2020162.aqua.250m.tiff"
falsecolor_test_image_file = "$(test_data_dir)/beaufort-chukchi-seas_falsecolor.2020162.aqua.250m.tiff"
falsecolor_b7_test_file = "$(test_data_dir)/ref_image_b7.png"
landmask_file = "$(test_data_dir)/landmask.tiff"
landmask_no_dilate_file = "$(test_data_dir)/julia_landmask_non_dilated.png"
current_landmask_file = "$(test_data_dir)/julia_landmask_dilated.png" # rename to julia dilated?
normalized_test_file = "$(test_data_dir)/matlab_normalized.png"
clouds_channel_test_file = "$(test_data_dir)/clouds_channel.png"
cloudmask_test_file = "$(test_data_dir)/cloudmask.png"
ice_water_discrim_test_file = "$(test_data_dir)/matlab_ice_water_discrim.png"
sharpened_test_image_file = "$(test_data_dir)/matlab_sharpened.png"
segmented_a_ice_mask_file = "$(test_data_dir)/matlab_segmented_A.png"
segmented_b_ice_test_file = "$(test_data_dir)/matlab_segmented_b_ice.png"
segmented_b_filled_test_file = "$(test_data_dir)/matlab_segmented_b_filled.png"
segmented_c_test_file = "$(test_data_dir)/matlab_segmented_c.png"
not_ice_mask_test_file = "$(test_data_dir)/matlab_not_ice_mask.png"
strel_file_2 = "$(test_data_dir)/se2.csv" # original matlab structuring element -  a disk-shaped kernel with radius of 2 px
strel_file_4 = "$(test_data_dir)/strel_disk_4.csv" # disk-shaped kernel with radius of 4 px
watershed_test_file = "$(test_data_dir)/matlab_watershed_intersect.png"

test_region = (1:2707, 1:4458)
lm_test_region = (1:800, 1:1500)
ice_floe_test_region = (1640:2060, 1840:2315)
