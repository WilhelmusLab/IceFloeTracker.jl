@testset "Discriminate Ice-Water" begin
    println("------------------------------------------------")
    println("------------ Create Discrimination Test --------------")

    Z, Z2 = IceFloeTracker.discriminate_ice_water(normalized_image)

    @test (@test_approx_eq_sigma_eps Z matlab_Z [0, 0] 0.005) == nothing

    @test (@test_approx_eq_sigma_eps Z2 matlab_Z2 [0, 0] 0.005) == nothing
end
