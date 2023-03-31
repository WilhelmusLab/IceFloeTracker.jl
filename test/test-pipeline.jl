@testset verbose = true "pipeline" begin
    using IceFloeTracker.Pipeline: cache_vector, extractfeatures, load_imgs, sharpen,sharpen_gray, disc_ice_water, load_truecolor_imgs, load_reflectance_imgs, cloudmask, get_ice_labels
    println("-------------------------------------------------")
    println("------------ pipeline funcs tests ---------------")

    pipelinedir = joinpath(test_data_dir, "pipeline")

    # input dir
    input = joinpath(pipelinedir, "input_pipeline")

    # output dir
    output = mkpath(joinpath(pipelinedir, "output"))

    args_to_pass = Dict{Symbol,AbstractString}(zip([:input, :output], [input, output]))

    reflectance_images = load_imgs(; input=input, image_type=:reflectance)

    truecolor_images = load_imgs(; input=input, image_type=:truecolor)

    lm_expected =
        IceFloeTracker.Gray.(
            load(joinpath(pipelinedir, "expected", "generated_landmask.png"))
        ) .> 0

    lm_raw = load(joinpath(input, "landmask.tiff"))

    landmask_no_dilate = (IceFloeTracker.Gray.(lm_raw) .> 0)
    sharpened_imgs = sharpen(truecolor_images, landmask_no_dilate)
    sharpenedgray_imgs = sharpen_gray(sharpened_imgs, lm_expected)

    @testset verbose = true "preprocessing" begin
        @testset "landmask" begin
            println("-------------------------------------------------")
            println("------------ landmask creation tests ---------------")

            IceFloeTracker.Pipeline.landmask(; args_to_pass...)
            @test isfile(joinpath(output, "generated_landmask.jls"))
            
            # deserialize landmask
            landmasks = deserialize(joinpath(output, "generated_landmask.jls"))
            
            @test lm_expected == landmasks.dilated
            @test .!(IceFloeTracker.Gray.(lm_raw) .> 0) == landmasks.non_dilated

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
            @test cldmasks_expected == cloudmask(; args_to_pass...)
        end

        @testset "load images" begin
            @test length(reflectance_images) == 2

            @test length(truecolor_images) == 2

            @test all(size.(reflectance_images) .== size.(truecolor_images))

            @test load_reflectance_imgs(; input=input) == reflectance_images
            @test load_truecolor_imgs(; input=input) == truecolor_images
        end

        @testset "ice water discrimination" begin
            cloudmasks = map(create_cloudmask, reflectance_images)
            normalized_images = [IceFloeTracker.normalize_image(sharpened_img, sharpened_gray_img, lm_expected) for (sharpened_img, sharpened_gray_img) in zip(sharpened_imgs, sharpenedgray_imgs)]
            ice_water_discrim_imgs = disc_ice_water(
                reflectance_images, normalized_images, cloudmasks, lm_expected
            )
            @test length(ice_water_discrim_imgs) == 2
        end

        @testset "ice labels" begin
            ice_labels = get_ice_labels(reflectance_images, lm_expected)
            @test length(ice_labels) == 2
        end

        @testset "sharpen_gray" begin
            @test eltype(sharpenedgray_imgs[1]) == Gray{Float64}
            @test length(sharpenedgray_imgs) == 2
        end
    end

    include("_test-preprocess.jl")

    @testset "feature extraction" begin
        minarea = "1"
        maxarea = "5"
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
                [:input, :output, :minarea, :maxarea, :features],
                [input, output, minarea, maxarea, features],
            ),
        )

        # generate two random image files with boolean data type using a seed
        container_to_serialize = cache_vector(Matrix{Bool}, 2, (200, 100))
        for i in eachindex(container_to_serialize)
            Random.seed!(i)
            container_to_serialize[i] .= .!rand((false, false, true, true, true), 200, 100)
        end

        serialize(joinpath(input, "segmented_floes.jls"), container_to_serialize)

        # run feature extraction
        @time extractfeatures(; args...)

        # check that the output files exist
        @test isfile(joinpath(output, "floe_props.jls"))

        # load the serialized output file
        floe_props = IceFloeTracker.deserialize(joinpath(output, "floe_props.jls"))
        @test typeof(floe_props) == Vector{DataFrame}
        @test length(floe_props) == 2

        # clean up!
        rm(extraction_path; recursive=true)
    end
end
