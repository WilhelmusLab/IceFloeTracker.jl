
@ntestset "$(@__FILE__)" begin
    @ntestset "Lopez-Acosta 2019" begin
        loadimg(s::String) = float64.(load(s))

        truecolor = loadimg(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        falsecolor = loadimg(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        landmask = loadimg("./test_inputs/pipeline/input_pipeline/landmask.tiff")

        results = LopezAcosta2019()(truecolor, falsecolor, landmask)
        # TODO: Add a check of the results
    end
    @ntestset "Lopez-Acosta 2019 with Tiling" begin
        truecolor = load(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        falsecolor = load(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        landmask = load("./test_inputs/pipeline/input_pipeline/landmask.tiff")

        results = LopezAcosta2019Tiling()(truecolor, falsecolor, landmask)
        # TODO: Add a check of the results
    end
end
