@testset verbose = true "pipeline" begin
    println("-------------------------------------------------")
    println("------------ pipeline funcs tests ---------------")

    pipelinedir = joinpath(test_data_dir, "pipeline")

    # input dir
    input = joinpath(pipelinedir, "input_pipeline")

    # output dir
    output = joinpath(pipelinedir, "output")

    args_to_pass = Dict{Symbol,AbstractString}(zip([:input, :output], [input, output]))

    @testset verbose = true "preprocessing" begin
        @testset "landmask" begin
            println("-------------------------------------------------")
            println("------------ landmask creation tests ---------------")
            lm_expected =
                Gray.(load(joinpath(pipelinedir, "expected", "generated_landmask.png"))) .>
                0

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
            reflectance_images = IceFloeTracker.load_imgs(;
                input=input, image_type=:reflectance
            )
            @test length(reflectance_images) == 2

            truecolor_images = IceFloeTracker.load_imgs(;
                input=input, image_type=:truecolor
            )
            @test length(truecolor_images) == 2

            @test all(size.(reflectance_images) .== size.(truecolor_images))

            @test IceFloeTracker.load_reflectance_imgs(; input=input) == reflectance_images
            @test IceFloeTracker.load_truecolor_imgs(; input=input) == truecolor_images
        end
    end
    # clean up!
    rm(output; recursive=true)
end
