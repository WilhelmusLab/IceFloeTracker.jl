@testset "Segmentation-D,E" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-D,E Test --------------")

    matlab_segmented_D = convert(
        BitMatrix,
        load("$(test_data_dir)/matlab_segmented_b_filled.png"),#
    )
    matlab_segmented_E = convert(
        BitMatrix,
        load("$(test_data_dir)/matlab_segmented_c.png"),#
    )
    matlab_watershed_D = convert(
        BitMatrix,
        load("$(test_data_dir)/matlab_watershed_B.png"),#
    )
    matlab_watershed_E = convert(
        BitMatrix,
        load("$(test_data_dir)/matlab_watershed_C.png"),#
    )
    matlab_watershed_intersect = convert(
        BitMatrix, load("$(test_data_dir)/matlab_watershed_intersect.png")
    )
    segmented_c = convert(BitMatrix, load(segmented_c_test_file))
    not_ice_mask = convert(BitMatrix, load(not_ice_mask_test_file))
    ## Run function with Matlab inputs
    @time watershed_D_borders = IceFloeTracker.segmentation_D(matlab_segmented_D) #not_ice_mask
    @time watershed_E_borders = IceFloeTracker.segmentation_E(matlab_segmented_E) #segmented_c
    @time watershed_intersect = IceFloeTracker.segmentation_D_E(
        watershed_D_borders, watershed_E_borders
    )

    watershed_D_filename =
        "$(test_output_dir)/watershed_d_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist watershed_D_borders watershed_D_filename

    watershed_E_filename =
        "$(test_output_dir)/watershed_e_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist watershed_E_borders watershed_E_filename

    watershed_intersect_filename =
        "$(test_output_dir)/watershed_intersect_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist watershed_intersect watershed_intersect_filename

    ## Run function with Julia inputs
    @time julia_watershed_D_borders = IceFloeTracker.segmentation_D(not_ice_mask) #not_ice_mask
    @time julia_watershed_E_borders = IceFloeTracker.segmentation_E(segmented_c) #segmented_c
    @time julia_watershed_intersect = IceFloeTracker.segmentation_D_E(
        julia_watershed_D_borders, julia_watershed_E_borders
    )

    julia_watershed_D_filename =
        "$(test_output_dir)/watershed_d_julia_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist julia_watershed_D_borders julia_watershed_D_filename

    julia_watershed_E_filename =
        "$(test_output_dir)/watershed_e_julia_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist julia_watershed_E_borders julia_watershed_E_filename

    julia_watershed_intersect_filename =
        "$(test_output_dir)/watershed_intersect_julia_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist julia_watershed_intersect julia_watershed_intersect_filename

    ## Tests with Matlab inputs
    @test typeof(watershed_E_borders) == typeof(matlab_watershed_D)
    @test typeof(watershed_D_borders) == typeof(matlab_watershed_E)
    @test typeof(watershed_intersect) == typeof(matlab_watershed_intersect)
    @test test_similarity(matlab_watershed_D, watershed_D_borders, 0.06)
    @test test_similarity(matlab_watershed_E, watershed_E_borders, 0.08)
    @test test_similarity(matlab_watershed_intersect, watershed_intersect, 0.017)

    ## Tests with Julia inputs
    @test typeof(julia_watershed_intersect) == typeof(matlab_watershed_intersect)
    @test test_similarity(matlab_watershed_D, julia_watershed_D_borders, 0.044)
    @test test_similarity(matlab_watershed_E, julia_watershed_E_borders, 0.09)
    @test test_similarity(matlab_watershed_intersect, julia_watershed_intersect, 0.009)
end
