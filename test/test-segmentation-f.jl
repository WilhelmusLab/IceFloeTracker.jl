@testset "Segmentation-F" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-F Test --------------")

    ## Load inputs for comparison
    segmentation_B_not_ice_mask = float64.(load("$(test_data_dir)/matlab_I.png"))
    matlab_Iobrd2 = float64.(load("$(test_data_dir)/matlab_Iobrd2.png"))
    matlab_Iobrcbr2 = float64.(load("$(test_data_dir)/matlab_Iobrcbr2.png"))
    matlab_C3_watershed = load("$(test_data_dir)/matlab_C3_watershed.png") .> 0.499
    matlab_BW1 = load("$(test_data_dir)/matlab_BW1.png") .> 0.499
    matlab_BW_final = load("$(test_data_dir)/matlab_BW_final.png") .> 0.499
    matlab_BW_final1 = load("$(test_data_dir)/matlab_BW_final1.png") .> 0.499
    matlab_BW_final2 = load("$(test_data_dir)/matlab_BW_final2.png") .> 0.499
    matlab_BW_final3 = load("$(test_data_dir)/matlab_BW_final3.png") .> 0.499
    matlab_BW_final4 = load("$(test_data_dir)/matlab_BW_final4.png") .> 0.499
    matlab_C4 = load("$(test_data_dir)/matlab_C4.png") .> 0.499
    matlab_C3_C4 = load("$(test_data_dir)/matlab_C3_C4.png") .> 0.499
    matlab_BW1_open = load("$(test_data_dir)/matlab_BW1_open.png") .> 0.499
    matlab_BW2 = load("$(test_data_dir)/matlab_BW2.png") .> 0.499
    matlab_BW4 = load("$(test_data_dir)/matlab_BW4.png") .> 0.499
    matlab_BW5 = load("$(test_data_dir)/matlab_BW5.png") .> 0.499
    matlab_BW6 = load("$(test_data_dir)/matlab_BW6.png") .> 0.499
    matlab_BW7 = load("$(test_data_dir)/matlab_BW7.png") .> 0.499
    segmentation_C_ice_mask = load(segmented_c_test_file) .> 0.499
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    landmask = convert(BitMatrix, load(current_landmask_file))
    watershed_intersect = load(watershed_test_file) .> 0.499
    ice_labels =
        Int64.(
            vec(DelimitedFiles.readdlm("$(test_data_dir)/ice_labels_floe_region.csv", ','))
        )
    
    ## Run function with Matlab inputs

    @time ice_mask_watershed_applied, ice_mask_watershed_opened, not_ice_dilated, reconstructed_leads, leads_segmented, leads_segmented_broken, 
    leads_branched, leads_filled, leads_opened_branched, leads_bothat, leads, leads_bothat_opened, leads_bothat_filled, leads_masked_branched, floes_erode, floes_dilate, floes_opened, isolated_floes = IceFloeTracker.segmentation_F(
        segmentation_C_ice_mask[ice_floe_test_region...],
        segmentation_B_not_ice_mask[ice_floe_test_region...],
        watershed_intersect[ice_floe_test_region...],
        cloudmask[ice_floe_test_region...],
        landmask[ice_floe_test_region...],
        ice_labels
    )
  
    IceFloeTracker.@persist isolated_floes "./test_outputs/isolated_floes.png" true
    
    IceFloeTracker.@persist not_ice_dilated "./test_outputs/not_ice_dilated.png" true
    IceFloeTracker.@persist matlab_Iobrd2[ice_floe_test_region...] "./test_outputs/matlab_Iobrd2.png" true

    IceFloeTracker.@persist (reconstructed_leads .> 0.499) "./test_outputs/reconstructed_leads.png" true

    IceFloeTracker.@persist matlab_Iobrcbr2[ice_floe_test_region...] .> 0.499 "./test_outputs/matlab_Iobrcbr2.png" true

    @test test_similarity(matlab_C3_watershed[ice_floe_test_region...], ice_mask_watershed_applied, 0.005)

    @test test_similarity(matlab_BW1[ice_floe_test_region...], ice_mask_watershed_opened, 0.006)

    @test (@test_approx_eq_sigma_eps matlab_Iobrd2[ice_floe_test_region...] not_ice_dilated [0, 0] 0.08) == nothing

    @test (@test_approx_eq_sigma_eps matlab_Iobrcbr2[ice_floe_test_region...] reconstructed_leads [0, 0] 0.135) == nothing

    @test test_similarity(matlab_BW_final[ice_floe_test_region...], leads_segmented, 0.035)

    @test test_similarity(matlab_BW_final1[ice_floe_test_region...], leads_segmented_broken, 0.035)

    @test test_similarity(matlab_BW_final2[ice_floe_test_region...], leads_branched, 0.047)

    @test test_similarity(matlab_BW_final3[ice_floe_test_region...], leads_filled, 0.047)

    @test test_similarity(matlab_BW_final4[ice_floe_test_region...], leads_opened_branched, 0.043)

    @test test_similarity(matlab_C4[ice_floe_test_region...], leads_bothat, 0.017)

    @test test_similarity(matlab_C3_C4[ice_floe_test_region...], leads, 0.043)

    @test test_similarity(matlab_BW1_open[ice_floe_test_region...], leads_bothat_opened, 0.033)

    @test test_similarity(matlab_BW2[ice_floe_test_region...], leads_bothat_filled, 0.068)

    IceFloeTracker.@persist leads_bothat_filled "./test_outputs/leads_bothat_filled.png" true
    
    IceFloeTracker.@persist matlab_BW2[ice_floe_test_region...] "./test_outputs/matlab_BW2_floe_region.png" true

    @test test_similarity(matlab_BW4[ice_floe_test_region...], leads_masked_branched, 0.05)

    @test test_similarity(matlab_BW5[ice_floe_test_region...], floes_erode, 0.027)

    @test test_similarity(matlab_BW6[ice_floe_test_region...], floes_opened, 0.052)

    @test test_similarity(matlab_BW7[ice_floe_test_region...], isolated_floes, 0.052)
end
