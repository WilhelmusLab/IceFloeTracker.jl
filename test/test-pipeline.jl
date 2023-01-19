@testset verbose = true "pipeline" begin
    println("-------------------------------------------------")
    println("------------ pipeline funcs tests ---------------")

    pipelinedir = joinpath(test_data_dir, "pipeline")

    # input dir
    input = joinpath(pipelinedir, "input_pipeline")

    # output dir
    output = joinpath(pipelinedir, "output")

    args_to_pass = Dict{Symbol,AbstractString}(zip([:input, :output], [input, output]))

    reflectance_images = IceFloeTracker.load_imgs(; input=input, image_type=:reflectance)

    truecolor_images = IceFloeTracker.load_imgs(; input=input, image_type=:truecolor)

    lm_expected =
        Gray.(load(joinpath(pipelinedir, "expected", "generated_landmask.png"))) .> 0

    @testset verbose = true "preprocessing" begin
        @testset "landmask" begin
            println("-------------------------------------------------")
            println("------------ landmask creation tests ---------------")

            @test lm_expected == IceFloeTracker.landmask(; args_to_pass...)

            @test isfile(joinpath(output, "generated_landmask.png"))
        end

        @testset "cloudmask" begin
            println("-------------------------------------------------")
            println("------------ cloudmasking tests -----------------")
            # Generate cloudmasks by hand
            cldmasks_paths = [f for f in readdir(input) if contains(f, "reflectance")]
            cldmasks_expected =
                IceFloeTracker.create_cloudmask.([
                    float64.(load(joinpath(input, f))) for f in cldmasks_paths
                ])

            # Compare against generated cloudmasks
            @test cldmasks_expected == IceFloeTracker.cloudmask(; args_to_pass...)
        end

        @testset "load images" begin
            @test length(reflectance_images) == 2

            @test length(truecolor_images) == 2

            @test all(size.(reflectance_images) .== size.(truecolor_images))

            @test IceFloeTracker.load_reflectance_imgs(; input=input) == reflectance_images
            @test IceFloeTracker.load_truecolor_imgs(; input=input) == truecolor_images
        end

        @testset "ice water discrimination" begin
            landmask_raw = load(joinpath(input, "landmask.tiff"))
            landmask_no_dilate = (Gray.(landmask_raw) .> 0)
            sharpened_imgs = IceFloeTracker.sharpen(truecolor_images, landmask_no_dilate)
            cloudmasks = map(create_cloudmask, reflectance_images)
            ice_water_discrim_imgs = IceFloeTracker.disc_ice_water(
                reflectance_images, sharpened_imgs, cloudmasks, lm_expected
            )
            @test length(ice_water_discrim_imgs) == 2
        end
    end
    # clean up!
    rm(output; recursive=true)
end
