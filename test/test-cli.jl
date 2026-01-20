@testitem "Segmentation algorithm filename call" begin
    using IceFloeTracker
    using FileIO
    using Images
    tempdir = mktempdir()
    # Download the dataset (default ref, cache_dir)
    dataset = Watkins2026Dataset()
    # Filter for cases with visible floes, no clouds, no artifacts
    filtered = filter(
        c -> (
            c.visible_floes == "yes" &&
            c.cloud_category_manual == "none" &&
            c.artifacts == "no"
        ),
        dataset,
    )
    # Take the first case for testing
    case = first(filtered)
    # Load required images
    truecolor = modis_truecolor(case)
    falsecolor = modis_falsecolor(case)
    landmask = modis_landmask(case)
    # Save images to temp files
    truecolor_path = joinpath(tempdir, "truecolor.tiff")
    falsecolor_path = joinpath(tempdir, "falsecolor.tiff")
    landmask_path = joinpath(tempdir, "landmask.tiff")
    output_path = joinpath(tempdir, "output.tiff")
    save(truecolor_path, truecolor)
    save(falsecolor_path, falsecolor)
    save(landmask_path, landmask)
    @testset "Run segmentation using file paths" begin
        seg = LopezAcosta2019Tiling.Segment()
        segmented = seg(output_path, truecolor_path, falsecolor_path, landmask_path)
        # Check output type and shape
        @test typeof(segmented) <: Images.SegmentedImage
        @test size(labels_map(segmented)) == size(truecolor)
    end

    @testset "Include an intermediate results callback" begin
        seg = LopezAcosta2019Tiling.Segment()
        segmented = seg(
            output_path,
            truecolor_path,
            falsecolor_path,
            landmask_path;
            intermediates_directory=tempdir,
            intermediates_targets=["icemask.tiff"],
        )
        # Check output type and shape
        @test typeof(segmented) <: Images.SegmentedImage
        @test size(labels_map(segmented)) == size(truecolor)
        @test isfile(joinpath(tempdir, "icemask.tiff"))
    end
end
