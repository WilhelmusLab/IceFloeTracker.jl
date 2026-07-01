@testitem "unknown formats aren't loaded" begin
    using NCDatasets
    mktemp() do output_path, _
        NCDataset(output_path, "c") do ds
            ds.attrib["file_version"] = "0.0.0"
            return nothing
        end
        @test_throws "file version" Archive.load(output_path)
    end
end

@testitem "choose_dtype" begin
    using IceFloeTracker.Archive: choose_dtype

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

@testsnippet ArchiveV1 begin
    using TimeZones
    using IceFloeTracker.Data:
        Watkins2026Dataset,
        modis_truecolor,
        modis_falsecolor,
        modis_landmask,
        modis_truecolor_path,
        modis_cloudfraction,
        pass_time,
        validated_labeled_floes,
        validated_floe_properties
    using FileIO
    using Images
    using DataFrames
    using NCDatasets

    dataset = Watkins2026Dataset()
    case = first(dataset)
    mask = modis_landmask(case)
    props = select(validated_floe_properties(case), Not(:boundary))
    # add the stereographic coordinate columns that the real pipeline produces
    props[!, :x] = zeros(nrow(props))
    props[!, :y] = zeros(nrow(props))
    data = IceFloeTracker.Archive.V1(;
        passtime=ZonedDateTime(pass_time(case), tz"UTC"),
        crs=latlon(modis_truecolor_path(case)),
        modis_truecolor=modis_truecolor(case),
        modis_falsecolor=modis_falsecolor(case),
        modis_cloud=RGB.(modis_cloudfraction(case)),
        labeled=validated_labeled_floes(case) |> labels_map,
        props=props,
        landmask=mask,
        cloud_mask=mask,
        ice_mask=mask,
        coastal_buffer_mask=mask,
        ift_version=VersionNumber("0.0.0"),
        reference="https://doi.org/00.0000",
        contact="contact@example.com",
        creation_date=ZonedDateTime(2024, 6, 15, 12, 6, 3, tz"UTC"),
        ift_configuration="IceFloeTracker.jl v0.0.0\n\nLopezAcosta2019Tiling.Segment\n",
    )
end

@testitem "Archive.V1 saved has the right global attributes" setup = [ArchiveV1] begin
    mktemp() do output_path, _
        Archive.save(output_path, data)
        NCDataset(output_path, "r") do ds
            @test ds.attrib["file_version"] == "1.0.0"
            @test ds.attrib["ift_version"] == "0.0.0"
            @test ds.attrib["reference"] == "https://doi.org/00.0000"
            @test ds.attrib["contact"] == "contact@example.com"
            @test ds.attrib["creation_date"] == "2024-06-15T12:06:03+00:00"
            @test ds.attrib["ift_configuration"] ==
                "IceFloeTracker.jl v0.0.0\n\nLopezAcosta2019Tiling.Segment\n"
        end
    end
end

@testitem "Archive.V1 saved has the right variables" setup = [ArchiveV1] begin
    mktemp() do output_path, _
        Archive.save(output_path, data)
        NCDataset(output_path, "r") do ds
            # coordinate variables
            @test haskey(ds, "x")
            @test haskey(ds, "y")
            @test haskey(ds, "time")
            @test haskey(ds, "geolocation")
            # colour imagery
            @test haskey(ds, "modis_truecolor")
            @test haskey(ds, "modis_falsecolor")
            @test haskey(ds, "modis_cloud")
            # segmentation outputs
            @test haskey(ds, "labeled_image")
            @test haskey(ds, "cloud_mask")
            @test haskey(ds, "landmask")
            @test haskey(ds, "ice_mask")
            @test haskey(ds, "coastal_buffer_mask")
            # floe properties
            @test haskey(ds, "floe_label")
            @test haskey(ds, "floe_area")
            @test haskey(ds, "floe_convex_area")
            # x/y prop columns must be remapped to avoid clashing with dimension names
            @test haskey(ds, "floe_x")
            @test haskey(ds, "floe_y")
            @test !haskey(ds, "x_crs")
        end
    end
end

@testitem "Archive.V1 truecolor and falsecolor have separate band dimensions" setup = [
    ArchiveV1
] begin
    mktemp() do output_path, _
        Archive.save(output_path, data)
        NCDataset(output_path, "r") do ds
            @test "band_modis_truecolor" in dimnames(ds["modis_truecolor"])
            @test "band_modis_falsecolor" in dimnames(ds["modis_falsecolor"])
            @test !("band_modis_falsecolor" in dimnames(ds["modis_truecolor"]))
            @test !("band_modis_truecolor" in dimnames(ds["modis_falsecolor"]))
        end
    end
end

@testitem "Archive.V1 can be saved with NaNs in the props" setup = [ArchiveV1] begin
    data.props[2, :convex_area] = NaN
    mktemp() do output_path, _
        @test_nowarn Archive.save(output_path, data)
    end
end

@testitem "Archive.V1 can be saved with missing in DataFrame, which converts to NaN" setup = [
    ArchiveV1
] begin
    allowmissing!(data.props, :convex_area)
    data.props[2, :convex_area] = missing
    mktemp() do output_path, _
        @test_nowarn Archive.save(output_path, data)
    end
end

@testitem "Archive.V1 saved can be reloaded correctly" setup = [ArchiveV1] begin
    mktemp() do output_path, _
        Archive.save(output_path, data)
        reloaded = Archive.load(output_path)
        @test reloaded.passtime == data.passtime
        @test reloaded.crs[:crs_wkt] == data.crs[:crs_wkt]
        @test reloaded.crs[:crs] == data.crs[:crs]
        @test reloaded.ift_version == data.ift_version
        @test reloaded.reference == data.reference
        @test reloaded.contact == data.contact
        @test reloaded.creation_date == data.creation_date
        @test reloaded.modis_truecolor == data.modis_truecolor
        @test reloaded.modis_falsecolor == data.modis_falsecolor
        @test reloaded.modis_cloud == data.modis_cloud
        @test reloaded.labeled == data.labeled
        @test reloaded.cloud_mask == Bool.(data.cloud_mask)
        @test reloaded.landmask == Bool.(data.landmask)
        @test reloaded.ice_mask == Bool.(data.ice_mask)
        @test reloaded.coastal_buffer_mask == Bool.(data.coastal_buffer_mask)
        @test isequal(reloaded.props, data.props)
        @test reloaded.ift_configuration == data.ift_configuration
    end
end
