@testset "Segmentation-D,E" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-D,E Test --------------")

    matlab_segmented_B = load("$(test_data_dir)/matlab_not_ice_mask.png") .<= 1
    matlab_segmented_C = convert(BitMatrix, load("$(test_data_dir)/matlab_segmented_c.png"))
    matlab_watershed_D = convert(BitMatrix, load("$(test_data_dir)/matlab_watershed_D.png"))
    matlab_watershed_E = convert(BitMatrix, load("$(test_data_dir)/matlab_watershed_E.png"))
    matlab_watershed_intersect = convert(
        BitMatrix, load("$(test_data_dir)/matlab_watershed_intersect.png")
    )

    ## Run function with Matlab inputs
    @time watershed_D_borders = IceFloeTracker.segmentation_D(matlab_segmented_C) #Matlab_segmented_c
    @time watershed_E_borders = IceFloeTracker.segmentation_E(matlab_segmented_B) #Matlab_not_ice_mask
    @time watershed_intersect = IceFloeTracker.segmentation_D_E(
        watershed_D_borders, watershed_E_borders
    )

    IceFloeTracker.@persist watershed_D_borders "./test_outputs/watershed_d.png" true
    
    IceFloeTracker.@persist watershed_D_borders "./test_outputs/watershed_e.png" true

    IceFloeTracker.@persist watershed_D_borders "./test_outputs/watershed_intersect.png" true

    ## Tests with Matlab inputs
    @test typeof(watershed_E_borders) == typeof(matlab_watershed_E)
    @test typeof(watershed_D_borders) == typeof(matlab_watershed_D)
    @test typeof(watershed_intersect) == typeof(matlab_watershed_intersect)
    @test test_similarity(matlab_watershed_D, watershed_D_borders, 0.078)
    @test test_similarity(matlab_watershed_E, watershed_E_borders, 0.024)
    @test test_similarity(matlab_watershed_intersect, watershed_intersect, 0.0028)
end
