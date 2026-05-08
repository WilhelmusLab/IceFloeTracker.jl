
@testitem "make_hdf5" begin
    using TimeZones
    using IceFloeTracker.Data:
        Watkins2026Dataset,
        modis_truecolor,
        modis_falsecolor,
        modis_landmask,
        modis_truecolor_path,
        modis_falsecolor_path,
        modis_landmask_path,
        pass_time
    using FileIO
    using Images
    using DataFrames

    mktemp() do output_path, _
        @info "Getting input files from validation dataset"
        dataset = Watkins2026Dataset()
        case = first(dataset)

        @show validated_floe_properties(case)

        make_hdf5(
            output_path;
            passtime=ZonedDateTime(pass_time(case), tz"UTC"),
            crs_ref_image_path=modis_truecolor_path(case),
            truecolor_path=modis_truecolor_path(case),
            falsecolor_path=modis_falsecolor_path(case),
            labeled=validated_labeled_floes(case) |> labels_map,
            props=select(validated_floe_properties(case), Not(:boundary)),
            landmask=modis_landmask(case),
            cloud_mask=modis_landmask(case), # placeholder for another mask type
            ice_mask=modis_landmask(case), # placeholder for another mask type
            coastal_buffer_mask=modis_landmask(case), # placeholder for another mask type
        )

        @test isfile(output_path)
    end
end