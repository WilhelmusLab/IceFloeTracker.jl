@ntestset "$(@__FILE__)" begin
    @ntestset "Lopez-Acosta 2019" begin
        loadimg(s::String) = float64.(load(s))

        truecolor = loadimg(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.truecolor.250m.tiff"
        )

        falsecolor = loadimg(
            "./test_inputs/pipeline/input_pipeline/20220914.aqua.reflectance.250m.tiff"
        )

        landmask = loadimg("./test_inputs/pipeline/input_pipeline/landmask.tiff")

        results = LopezAcosta2019()(truecolor, falsecolor, landmask)
        @show results
    end
end
