@testsnippet HDF5V1 begin
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
    using HDF5

    dataset = Watkins2026Dataset()
    case = first(dataset)
    data = IceFloeTracker.HDF5.V1(;
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
        iftversion=VersionNumber("0.0.0"),
        reference="https://doi.org/00.0000",
        contact="contact@example.com",
    )
end

@testitem "HDF5.V1 saved have the right fields" setup = [HDF5V1] begin
    mktemp() do output_path, _
        save_hdf5(output_path, data;)
        h5open(output_path, "r") do file
            @test attrs(file)["file_version"] === "1.0.0"
            @test attrs(file)["iftversion"] === "0.0.0"
            @test attrs(file)["reference"] === "https://doi.org/00.0000"
            @test attrs(file)["contact"] === "contact@example.com"
            @test attrs(file)["fname_truecolor"] == modis_truecolor_path(case)
            @test attrs(file)["crs_name"] === "EPSG:3413 NSIDC north polar stereographic"
            @test haskey(file, "index")
            @test haskey(file, "floe_properties")
            @test haskey(file["floe_properties"], "labeled_image")
            @test haskey(file["floe_properties"], "properties")
            @test haskey(file, "classifications")
            @test haskey(file["classifications"], "landmask")
            @test haskey(file["classifications"], "ice_mask")
            @test haskey(file["classifications"], "coastal_buffer_mask")
        end
    end
end

@testitem "HDF5.V1 saved can be reloaded correctly" setup = [HDF5V1] begin
    mktemp() do output_path, _
        save_hdf5(output_path, data;)
        reloaded = load_hdf5(output_path)
        @test reloaded.passtime == data.passtime
        @test reloaded.crs_ref_image_path == data.crs_ref_image_path
        @test reloaded.truecolor_path == data.truecolor_path
        @test reloaded.falsecolor_path == data.falsecolor_path
        @test reloaded.labeled == data.labeled
        @test isequal(reloaded.props, data.props)
        @test reloaded.cloud_mask == data.cloud_mask
        @test reloaded.ice_mask == data.ice_mask
        @test reloaded.landmask == data.landmask
        @test reloaded.coastal_buffer_mask == data.coastal_buffer_mask
        @test reloaded.iftversion == data.iftversion
        @test reloaded.reference == data.reference
        @test reloaded.contact == data.contact
    end
end

@testitem "HDF5.V1 can be saved with NaNs in the props" setup = [HDF5V1] begin
    @show data.props
    data.props[2, :convex_area] = NaN
    @show data.props
    mktemp() do output_path, _
        save_hdf5(output_path, data;)
        reloaded = load_hdf5(output_path)
        @test isequal(reloaded.props, data.props)
    end
end

@testitem "HDF5.V1 can be saved with missing in DataFrame, which convert to NaN" setup = [
    HDF5V1
] begin
    convex_area = allowmissing!(data.props, :convex_area)
    data.props[2, :convex_area] = missing
    mktemp() do output_path, _
        save_hdf5(output_path, data;)
        reloaded = load_hdf5(output_path)
        @test isequal(reloaded.props[1, :], data.props[1, :])
        @test isequal(reloaded.props[2, :convex_area], NaN)  # missing becomes a NaN when saved and reloaded
        @test isequal(reloaded.props[3:end, :], data.props[3:end, :])
    end
end

@testitem "unknown HDF5 files aren't loaded" begin
    using HDF5

    mktemp() do output_path, _
        h5open(output_path, "w") do file
            return attrs(file)["file_version"] = "0.0.0"
        end
        @test_throws "file version" load_hdf5(output_path)
    end
end

@testitem "choose_dtype" begin
    using IceFloeTracker.HDF5: choose_dtype

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

@testitem "convert_missing_to_nan" begin
    using DataFrames
    using IceFloeTracker.HDF5: convert_missing_to_nan

    df = DataFrame(a=[1.0, 2.0, missing], b=[missing, 4.0, 5.0], c=[6.0, 7.0, 8.0])
    df_converted = convert_missing_to_nan(df)
    @test isequal(df_converted.a, [1.0, 2.0, NaN])
    @test isequal(df_converted.b, [NaN, 4.0, 5.0])
    @test isequal(df_converted.c, [6.0, 7.0, 8.0])
end