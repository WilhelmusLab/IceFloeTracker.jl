@testset "register_mismatch tests" begin
    # Read in floe from file
    floe = readdlm("./test_inputs/floetorotate.csv", ',', Bool)

    # Define rotation angle
    rot_angle = 15 # degrees

    # Rotate floe
    rotated_floe = imrotate(floe, deg2rad(rot_angle))

    # Estimate rigid transformation and mismatch
    mm, rot = IceFloeTracker.mismatch(floe, rotated_floe)

    # Test 1: mismatch accuracy
    @test mm < 0.005

    # Test 2: angle estimate
    @test abs(rot - rot_angle) < 0.05
end
