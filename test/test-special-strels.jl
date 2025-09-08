
@testitem "Special strels" begin
    using IceFloeTracker: se_disk4, se_disk20, se_disk50
    @testset "se_disk4" begin
        @test sum(se_disk4()) == 37
    end

    @testset "se_disk20" begin
        @test sum(se_disk20()) == 1301
    end

    @testset "se_disk50/se for landmask dilation" begin
        @test sum(se_disk50()) == 8177
    end
end
