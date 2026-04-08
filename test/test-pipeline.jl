@testsnippet Pipeline begin
    function pipeline_runs(pipeline_module, dataset::Dataset)::Bool
        truecolor = modis_truecolor.(dataset)
        falsecolor = modis_falsecolor.(dataset)
        landmask = modis_landmask.(dataset)

        segmentation = pipeline_module.Segment().(truecolor, falsecolor, landmask)
        pass_times = pass_time.(dataset)  # Assuming passtime is a function that extracts the time from each case

        tracking = pipeline_module.Track()(segmentation, pass_times)
        # Add any additional checks here if needed, e.g., ensuring that tracking results are not empty or have expected properties

        return true  # If we reach this point without errors, the test passes
    end
end
@testitem "Pipelines" tags = [:e2e] setup = [Pipeline] begin
    dataset = Watkins2026Dataset(; ref="v0.2")
    @test pipeline_runs(LopezAcosta2019, filter(c -> c.case_number == 4, dataset))
    @test pipeline_runs(LopezAcosta2019Tiling, filter(c -> c.case_number == 4, dataset))
end
