@testitem "Segmentation-Watershed" begin
    include("config.jl")
    include("test_error_rate.jl")

    matlab_not_ice = load("$(test_data_dir)/matlab_not_ice_mask.png")
    matlab_not_ice_bit = matlab_not_ice .> 0.499
    matlab_ice_intersect = convert(
        BitMatrix, load("$(test_data_dir)/matlab_segmented_c.png")
    )
    matlab_watershed_D = convert(BitMatrix, load("$(test_data_dir)/matlab_watershed_D.png"))
    matlab_watershed_E = convert(BitMatrix, load("$(test_data_dir)/matlab_watershed_E.png"))
    matlab_watershed_intersect = convert(
        BitMatrix, load("$(test_data_dir)/matlab_watershed_intersect.png")
    )

    ## Run function with Matlab inputs
    @time watershed_B_ice_intersect = IceFloeTracker.watershed_ice_floes(
        matlab_ice_intersect
    )
    @time watershed_B_not_ice = IceFloeTracker.watershed_ice_floes(matlab_not_ice_bit)
    @time watershed_intersect = IceFloeTracker.watershed_product(
        watershed_B_ice_intersect, watershed_B_not_ice
    )

    IceFloeTracker.@persist watershed_B_ice_intersect "./test_outputs/watershed_ice_intersect.png" true

    IceFloeTracker.@persist watershed_B_not_ice "./test_outputs/watershed_not_ice.png" true

    IceFloeTracker.@persist watershed_intersect "./test_outputs/watershed_intersect.png" true

    IceFloeTracker.@persist matlab_not_ice_bit "./test_outputs/matlab_not_ice_bit.png" true

    ## Tests with Matlab inputs
    @test typeof(watershed_B_not_ice) == typeof(matlab_watershed_E)
    @test typeof(watershed_B_ice_intersect) == typeof(matlab_watershed_D)
    @test typeof(watershed_intersect) == typeof(matlab_watershed_intersect)
    @test test_similarity(matlab_watershed_D, watershed_B_ice_intersect, 0.12)
    @test test_similarity(matlab_watershed_E, watershed_B_not_ice, 0.15)
    @test test_similarity(matlab_watershed_intersect, watershed_intersect, 0.033)
end
