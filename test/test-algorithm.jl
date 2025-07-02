loadimg(s::String) = (x -> float64.(x))(load(s))

@testset "algorithm tests" begin
    println("-------------------------------------------------")
    println("----------- high-level algorithm tests ----------")
    @testset "Lopez-Acosta 2019" begin
        truecolor = loadimg(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        falsecolor = loadimg(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        landmask = loadimg("./test_inputs/pipeline/input_pipeline/landmask.tiff")

        results = LopezAcosta2019()(truecolor, falsecolor, landmask)
        @info results
    end
    @testset "Lopez-Acosta 2019 with Tiling" begin
        truecolor = loadimg(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        falsecolor = loadimg(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        landmask = loadimg("./test_inputs/pipeline/input_pipeline/landmask.tiff")

        results = LopezAcosta2019Tiling()(truecolor, falsecolor, landmask)
        @info results
    end
end
