function test_cloud_image_workflow()
    @testset "Prereq cloud image" begin

        datadir = joinpath(@__DIR__, "test_inputs/")
        path_ref_image = joinpath(datadir, "NE_Greenland_reflectance.2020162.aqua.250m.tiff")
        ref_image = float64.(load(path_ref_image))
        dilated_landmask = BitMatrix(load(joinpath(datadir, "matlab_landmask.png")))

        redchannel_cahe = IceFloeTracker._get_red_channel_cloud_cae(
            false_color_image=ref_image,
            landmask=dilated_landmask,
            prelim_threshold=110.0,
            band_7_threshold=200.0,
            band_2_threshold=190.0,
        )

        @test sum(redchannel_cahe) == 1320925065
    end
end


function test_adaphisteq()
    @testset "Adaptive histogram equalization" begin
        img = IceFloeTracker.convert_to_255_matrix(testimage("cameraman"))
        img_eq = IceFloeTracker.adapthisteq(img)
        @test sum(img_eq) == 32387397
    end
end

@testset "Conditional adaptivehisteq" begin

    test_cloud_image_workflow()
    test_adaphisteq()

end
