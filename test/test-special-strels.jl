
@testitem "Special strels" begin
    import IceFloeTracker: strel_octagon, make_landmask_se, strel_disk
    import Images: strel_diamond

    @testset "Octagon with radius 3 (se_disk4)" begin
        @test strel_octagon(3) == centered(
            Bool[
                0 0 1 1 1 0 0
                0 1 1 1 1 1 0
                1 1 1 1 1 1 1
                1 1 1 1 1 1 1
                1 1 1 1 1 1 1
                0 1 1 1 1 1 0
                0 0 1 1 1 0 0
            ],
        )
    end

    @testset "Octagon with radius 19 (se_disk20)" begin
        @test sum(strel_octagon(19)) == 1301
    end

    @testset "Special case for landmask dilation" begin
        @test sum(make_landmask_se()) == 8177
    end

    @testset "Small disk is the same as a diamond" begin
        @test strel_disk(1) == strel_diamond((3,3))
    end

    @testset "Large disk is approximately circular" begin
        @test sum(strel_disk(400)) ≈ (pi*400^2) atol = 1e-3
    end


end
