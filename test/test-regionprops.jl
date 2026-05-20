@testitem "regionprops" begin
    using Random
    import DataFrames: DataFrame, nrow
    import Images: label_components
    import IceFloeTracker: PixelConvexArea, PolygonConvexArea

    Random.seed!(123)
    bw_img = Bool.(rand([0, 1], 5, 10))
    bw_img[end, 7] = 1
    label_img = label_components(bw_img, trues(3, 3))
    properties = [
        "label",
        "centroid",
        "area",
        "major_axis_length",
        "minor_axis_length",
        "bbox",
        "perimeter",
        "orientation",
    ]
    extra_props = nothing

    table = regionprops_table(label_img, bw_img; properties=properties, minimum_area=1)
    total_labels = maximum(label_img)

    # Tests for regionprops_table
    @test typeof(table) <: DataFrame # check correct data type
    @test nrow(table) == total_labels - 1 # One point has area 1, we exclude that point intentionally
    @test all(
        [
            "area",
            "min_row",
            "min_col",
            "max_row",
            "max_col",
            "row_centroid",
            "col_centroid",
            "label",
            "major_axis_length",
            "minor_axis_length",
            "orientation",
            "perimeter",
        ] ⊆ names(table),
    )
    # check all requested properties are present

    # Tests for regionprops
    regions = regionprops(label_img, bw_img)
    @test table.area[2] == regions[:area][2]

    # Check floe masks generation and correct cropping
    add_floemasks!(table, label_img)
    @test all(
        [
            length(unique(label_components(table.mask[i], trues(3, 3)))) for
            i in 1:nrow(table)
        ] .== [2, 2],
    )

    # Test that algorithm options at least run
    regionprops(label_img; properties=["perimeter"], perimeter_algorithm=BenkridCrookes())
    regionprops(
        label_img; properties=["convex_area"], convex_area_algorithm=PolygonConvexArea()
    )
    regionprops(
        label_img; properties=["convex_area"], convex_area_algorithm=PixelConvexArea()
    )
end

@testitem "regionprops table output should include all implied columns even if there are no rows" begin
    using DataFrames: nrow

    label_img = Int.(zeros(5, 5))
    @testset "defaults" begin
        table = regionprops_table(label_img)
        @test "min_row" ∈ names(table)
        @test "max_row" ∈ names(table)
        @test "min_col" ∈ names(table)
        @test "max_col" ∈ names(table)
        @test "row_centroid" ∈ names(table)
        @test "col_centroid" ∈ names(table)
        @test "label" ∈ names(table)
        @test "area" ∈ names(table)
        @test "major_axis_length" ∈ names(table)
        @test "minor_axis_length" ∈ names(table)
        @test "convex_area" ∈ names(table)
        @test "perimeter" ∈ names(table)
        @test "orientation" ∈ names(table)
    end

    @testset "bbox" begin
        table = regionprops_table(label_img; properties=["bbox"])
        @test nrow(table) == 0
        @test "min_row" ∈ names(table)
        @test "max_row" ∈ names(table)
        @test "min_col" ∈ names(table)
        @test "max_col" ∈ names(table)
    end

    @testset "centroid" begin
        table = regionprops_table(label_img; properties=["centroid"])
        @test "row_centroid" ∈ names(table)
        @test "col_centroid" ∈ names(table)
    end
end

@testitem "addlatlon! adds latitude and longitude columns to regionprops table" begin
    using IceFloeTracker
    import DataFrames: DataFrame
    using FileIO

    refimage = "test_inputs/latlon/latlon_test_image-2020-06-21T00_00_00Z.tif"
    data = load(refimage)
    nrows, ncols = size(data)

    propdf = DataFrame(;
        row_centroid=[1.0, 1.0, Float64(nrows), Float64(nrows)],
        col_centroid=[1.0, Float64(ncols), 1.0, Float64(ncols)],
        area=ones(4),
        convex_area=ones(4),
        minor_axis_length=ones(4),
        major_axis_length=ones(4),
        perimeter=ones(4),
    )

    addlatlon!(propdf, refimage)

    @test all(["latitude", "longitude"] ⊆ names(propdf)) # check new columns added
    @info propdf
    top_left, top_right, bottom_left, bottom_right = 1, 2, 3, 4

    # For this georeferenced image, northern pixels (top rows) should have higher latitude.
    @test propdf.latitude[top_left] > propdf.latitude[bottom_left]
    @test propdf.latitude[top_right] > propdf.latitude[bottom_right]

    # Longitudes should increase from left to right at both top and bottom rows.
    @test propdf.longitude[top_right] > propdf.longitude[top_left]
    @test propdf.longitude[bottom_right] > propdf.longitude[bottom_left]
end

@testitem "converttounits maps corner centroids to sensible lat/lon ordering" begin
    using IceFloeTracker
    import DataFrames: DataFrame

    refimage = "test_inputs/latlon/latlon_test_image-2020-06-21T00_00_00Z.tif"
    latlondata = latlon(refimage)
    latitude = latlondata[:latitude]
    nrows, ncols = size(latitude)

    propdf = DataFrame(;
        row_centroid=[1.0, 1.0, Float64(nrows), Float64(nrows)],
        col_centroid=[1.0, Float64(ncols), 1.0, Float64(ncols)],
        area=ones(4),
        convex_area=ones(4),
        minor_axis_length=ones(4),
        major_axis_length=ones(4),
        perimeter=ones(4),
    )

    converted = converttounits(propdf, latlondata)

    top_left, top_right, bottom_left, bottom_right = 1, 2, 3, 4

    # For this georeferenced image, northern pixels (top rows) should have higher latitude.
    @test converted.latitude[top_left] > converted.latitude[bottom_left]
    @test converted.latitude[top_right] > converted.latitude[bottom_right]

    # Longitudes should increase from left to right at both top and bottom rows.
    @test converted.longitude[top_right] > converted.longitude[top_left]
    @test converted.longitude[bottom_right] > converted.longitude[bottom_left]
end

@testitem "converttounits copes with empty data frames" begin
    using IceFloeTracker
    import DataFrames: DataFrame

    refimage = "test_inputs/latlon/latlon_test_image-2020-06-21T00_00_00Z.tif"
    latlondata = latlon(refimage)
    latitude = latlondata[:latitude]
    nrows, ncols = size(latitude)

    propdf = DataFrame(;
        row_centroid=[],
        col_centroid=[],
        area=Float64[],
        convex_area=Float64[],
        minor_axis_length=Float64[],
        major_axis_length=Float64[],
        perimeter=Float64[],
    )

    @test converttounits(propdf, latlondata) isa DataFrame # doesn't crash
end
