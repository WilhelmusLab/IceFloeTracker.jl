
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

@testitem "choose_dtype" begin
    using IceFloeTracker.PersistHDF5: choose_dtype

    @test choose_dtype(100) == UInt8
    @test choose_dtype(-100) == Int8
    @test choose_dtype(1000) == UInt16
    @test choose_dtype(-1000) == Int16
    @test choose_dtype(100000) == UInt32
    @test choose_dtype(-100000) == Int32
    @test choose_dtype(10000000000) == UInt64
    @test choose_dtype(-10000000000) == Int64
    @test choose_dtype(BigInt(2)^64 - 1) == UInt64
    @test choose_dtype(-BigInt(2)^63) == Int64
    @test_throws ErrorException choose_dtype(BigInt(2)^64 + 1)
    @test_throws ErrorException choose_dtype(-BigInt(2)^63 - 1)
end