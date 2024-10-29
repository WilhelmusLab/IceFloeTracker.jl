using IceFloeTracker: getlatlon

@testset "getlatlon" begin
    imgpth = "test_inputs/latlon/latlon_test_image-2020-06-21T00_00_00Z.tif"
    expected_lat = readdlm("test_inputs/latlon/latitude.csv", ',', Float64)
    expected_lon = readdlm("test_inputs/latlon/longitude.csv", ',', Float64)

    latlon = getlatlon(imgpth)
    @test expected_lat == latlon["latitude"]
    @test expected_lon == latlon["longitude"]
end
