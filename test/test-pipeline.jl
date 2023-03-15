@testset verbose = true "pipeline" begin
    println("-------------------------------------------------")
    println("------------ pipeline funcs tests ---------------")

    pipelinedir = joinpath(test_data_dir, "pipeline")

    # input dir
    input = joinpath(pipelinedir, "input_pipeline")

    # output dir
    output = mkpath(joinpath(pipelinedir, "output"))

    args_to_pass = Dict{Symbol,AbstractString}(zip([:input, :output], [input, output]))

    reflectance_images = IceFloeTracker.load_imgs(; input=input, image_type=:reflectance)

    truecolor_images = IceFloeTracker.load_imgs(; input=input, image_type=:truecolor)

    lm_expected =
        IceFloeTracker.Gray.(
            load(joinpath(pipelinedir, "expected", "generated_landmask.png"))
        ) .> 0

    lm_raw = load(joinpath(input, "landmask.tiff"))

    landmask_no_dilate = (IceFloeTracker.Gray.(lm_raw) .> 0)
    sharpened_imgs = IceFloeTracker.sharpen(truecolor_images, landmask_no_dilate)
    sharpenedgray_imgs = IceFloeTracker.sharpen_gray(sharpened_imgs, lm_expected)

    @testset verbose = true "preprocessing" begin
        @testset "landmask" begin
            println("-------------------------------------------------")
            println("------------ landmask creation tests ---------------")

            IceFloeTracker.landmask(; args_to_pass...)
            @test isfile(joinpath(output, "generated_landmask.jls"))
            
            # deserialize landmask
            landmasks = deserialize(joinpath(output, "generated_landmask.jls"))
            
            @test lm_expected == landmasks.dilated
            @test .!(IceFloeTracker.Gray.(lm_raw) .> 0) == landmasks.non_dilated
            @test isfile(joinpath(output, "generated_landmask_dilated.png"))
            @test isfile(joinpath(output, "generated_landmask_non_dilated.png"))

            # clean up!
            rm(output; recursive=true)
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
            cloudmasks = map(create_cloudmask, reflectance_images)
            normalized_images = [IceFloeTracker.normalize_image(sharpened_img, sharpened_gray_img, lm_expected) for (sharpened_img, sharpened_gray_img) in zip(sharpened_imgs, sharpenedgray_imgs)]
            ice_water_discrim_imgs = IceFloeTracker.disc_ice_water(
                reflectance_images, normalized_images, cloudmasks, lm_expected
            )
            @test length(ice_water_discrim_imgs) == 2
        end

        @testset "ice labels" begin
            ice_labels = IceFloeTracker.get_ice_labels(reflectance_images, lm_expected)
            @test length(ice_labels) == 2
        end

        @testset "sharpen_gray" begin
            @test eltype(sharpenedgray_imgs[1]) == Gray{Float64}
            @test length(sharpenedgray_imgs) == 2
        end
    end

    @testset "feature extraction" begin
        min_area = "1"
        max_area = "5"
        features = "area bbox centroid"
        extraction_path = joinpath(
            @__DIR__, "test_inputs", "pipeline", "feature_extraction"
        )
        ispath(joinpath(extraction_path, "input")) &&
            rm(joinpath(extraction_path, "input"); recursive=true)
        input = mkpath(joinpath(extraction_path, "input"))
        output = mkpath(joinpath(extraction_path, "output"))
        args = Dict{Symbol,Any}(
            zip(
                [:input, :output, :min_area, :max_area, :features],
                [input, output, min_area, max_area, features],
            ),
        )

        # generate two random image files with boolean data type using a seed
        for i in 1:2
            Random.seed!(i)
            @persist .!rand((false, false, true, true, true), 200, 100) joinpath(
                input, "floe$i.png"
            )
        end

        # run feature extraction
        IceFloeTracker.extractfeatures(; args...)

        # check that the output files exist
        @test isfile(joinpath(output, "floe_library.dat"))

        # load the serialized output file
        floe_library = IceFloeTracker.deserialize(joinpath(output, "floe_library.dat"))
        @test typeof(floe_library) == Vector{DataFrame}
        @test length(floe_library) == 2

        # clean up!
        rm(extraction_path; recursive=true)
    end
end
