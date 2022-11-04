@testset "register_mismatch tests" begin
    # Read in floe from file
    floe = readdlm("./test_inputs/floetorotate.csv", ',', Bool)

    # Define rotation angle
    rot_angle = deg2rad(15)

    # Rotate floe
    rotated_floe = imrotate(floe, rot_angle)

    # Estimate rigid transformation and mismatch
    mm, tfm = IceFloeTracker.mismatch(floe, rotated_floe)

    # Test 1: mismatch accuracy
    @test mm < 0.005

    # Test 2: angle estimate
    @test abs(acos(tfm.linear[1]) - rot_angle) < 0.05
end
