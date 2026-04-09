@testitem "Segmentation algorithm filename call" begin
    using IceFloeTracker
    using FileIO
    using Images
    using IceFloeTracker.Preprocessing: binarize_mask
    using IceFloeTracker.Utils: call_kwargs

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
    ice_mask_path = joinpath(tempdir, "ice_mask.tiff")
    save(truecolor_path, truecolor)
    save(falsecolor_path, falsecolor)
    save(landmask_path, landmask)

    function test_algorithm_with_callback(seg)
        segmented = seg(
            load(truecolor_path),
            load(falsecolor_path),
            load(landmask_path) |> binarize_mask;
            intermediate_results_callback=call_kwargs(;
                labels_map=l -> l .|> UInt16 |> save(output_path),
                ice_mask=save(ice_mask_path),
            ),
        )
        result = (;
            imageType=typeof(segmented) <: Images.SegmentedImage,
            imageSize=size(labels_map(segmented)) == size(truecolor),
            icemaskExists=isfile(ice_mask_path),
        )
        if !all(values(result))
            @warn "Intermediate results callback test failed" result = result
        end
        return result
    end

    @test all(test_algorithm_with_callback(LopezAcosta2019.Segment()))
    @test all(test_algorithm_with_callback(LopezAcosta2019Tiling.Segment()))
end
