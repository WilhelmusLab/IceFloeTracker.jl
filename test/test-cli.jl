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

    function test_algorithm(seg)
        segmented = seg(output_path, truecolor_path, falsecolor_path, landmask_path)
        # Check output type and shape
        return typeof(segmented) <: Images.SegmentedImage &&
               size(labels_map(segmented)) == size(truecolor)
    end

    @test test_algorithm(LopezAcosta2019.Segment())
    @test test_algorithm(LopezAcosta2019Tiling.Segment())

    function test_algorithm_with_callback(seg)
        segmented = seg(
            output_path,
            truecolor_path,
            falsecolor_path,
            landmask_path;
            intermediates_directory=tempdir,
            intermediates_targets=["icemask.tiff"],
        )
        result = (;
            imageType=typeof(segmented) <: Images.SegmentedImage,
            imageSize=size(labels_map(segmented)) == size(truecolor),
            icemaskExists=isfile(joinpath(tempdir, "icemask.tiff")),
        )
        if !all(values(result))
            @warn "Intermediate results callback test failed" result = result
        end
        return result
    end

    @test all(test_algorithm_with_callback(LopezAcosta2019.Segment()))
    @test all(test_algorithm_with_callback(LopezAcosta2019Tiling.Segment()))
end
