using IceFloeTracker: getlatlon

@testset "getlatlon" begin
    imgpth = "test_inputs/latlon/latlon_test_image-2020-06-21T00_00_00Z.tif"

    # vec needed to convert to vector instead of nx1 matrix
    expected_X = vec(readdlm("test_inputs/latlon/X.csv", ',', Float64))
    expected_Y = vec(readdlm("test_inputs/latlon/Y.csv", ',', Float64))

    expected_lat = readdlm("test_inputs/latlon/latitude.csv", ',', Float64)
    expected_lon = readdlm("test_inputs/latlon/longitude.csv", ',', Float64)

    latlon = getlatlon(imgpth)
    X = latlon["X"]
    Y = latlon["Y"]
    lat = latlon["latitude"]
    lon = latlon["longitude"]

    # round to 4 decimal places to avoid weirdness in different arch/os
    expected_lat, expected_lon, expected_X, expected_Y = [round.(arr, digits=4) for arr in [expected_lat, expected_lon, expected_X, expected_Y]]
    lat, lon, X, Y = [round.(arr, digits=4) for arr in [lat, lon, X, Y]]

    @test expected_X == X
    @test expected_Y == Y
    @test expected_lat == lat
    @test expected_lon == lon
end
