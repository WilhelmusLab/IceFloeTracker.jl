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
