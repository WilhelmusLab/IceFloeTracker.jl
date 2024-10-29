using IceFloeTracker: getlatlon

@testset "getlatlon" begin
    imgpth = "test_inputs/latlon/latlon_test_image-2020-06-21T00_00_00Z.tif"

    # vec needed to convert to vector instead of nx1 matrix
    expected_X = vec(readdlm("test_inputs/latlon/X.csv", ',', Float64))
    expected_Y = vec(readdlm("test_inputs/latlon/Y.csv", ',', Float64))

    expected_lat = readdlm("test_inputs/latlon/latitude.csv", ',', Float64)
    expected_lon = readdlm("test_inputs/latlon/longitude.csv", ',', Float64)
    
    latlon = getlatlon(imgpth)
    @test expected_X == latlon["X"]
    @test expected_Y == latlon["Y"]
    @test expected_lat == latlon["latitude"]
    @test expected_lon == latlon["longitude"]
end
