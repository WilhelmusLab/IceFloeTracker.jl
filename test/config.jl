# Setting things up

## locate some files for the tests
test_data_dir = "./test_inputs"
test_output_dir = "./test_outputs"
truecolor_test_image_file = "$(test_data_dir)/NE_Greenland_truecolor.2020162.aqua.250m.tiff"
reflectance_test_image_file = "$(test_data_dir)/NE_Greenland.2020162.aqua.250m.tiff"
reflectance_b7_test_file = "$(test_data_dir)/ref_image_b7.png"
landmask_file = "$(test_data_dir)/landmask.tiff"
landmask_no_dilate_file = "$(test_data_dir)/landmask_no_dilate.png"
current_landmask_file = "$(test_data_dir)/current_landmask.png"
normalized_test_file = "$(test_data_dir)/normalized_image.png"
clouds_channel_test_file = "$(test_data_dir)/clouds_channel.png"
cloudmask_test_file = "$(test_data_dir)/cloudmask.png"
ice_water_discrim_test_file = "$(test_data_dir)/ice_water_discrim_image.png"
sharpened_test_image_file = "$(test_data_dir)/sharpened_test_image.png"
segmented_a_ice_mask_file = "$(test_data_dir)/segmented_ice_cloudmasked.png"
segmented_b_ice_test_file = "$(test_data_dir)/segmented_B_ice.png"
segmented_b_filled_test_file = "$(test_data_dir)/segmented_B_filled.png"
segmented_c_test_file = "$(test_data_dir)/segmented_C.png"
not_ice_mask_test_file = "$(test_data_dir)/not_ice_mask.png"
strel_file_2 = "$(test_data_dir)/se2.csv" # original matlab structuring element -  a disk-shaped kernel with radius of 2 px

test_region = (1:2707, 1:4458)
lm_test_region = (1:800, 1:1500)
ice_floe_test_region = (1640:2060, 1840:2315)
