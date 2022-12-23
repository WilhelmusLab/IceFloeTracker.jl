@testset verbose = true "pipeline" begin
    println("-------------------------------------------------")
    println("------------ pipeline funcs tests ---------------")

    pipelinedir = joinpath(test_data_dir, "pipeline")

    # input dir
    input = joinpath(pipelinedir, "input_pipeline")

    # output dir
    output = joinpath(pipelinedir, "output")

    @testset verbose = true "preprocessing" begin
        @testset "landmask" begin
            println("-------------------------------------------------")
            println("------------ landmask creation tests ---------------")
            lm_expected =
                Gray.(load(joinpath(pipelinedir, "expected", "generated_landmask.png"))) .>
                0
            args_to_pass = Dict{Symbol,AbstractString}(
                zip([:input, :output], [input, output])
            )
            @test lm_expected == IceFloeTracker.landmask(; args_to_pass...)

            @test isfile(joinpath(output, "generated_landmask.png"))
        end
    end

    # clean up!
    rm(output; recursive=true)
end
