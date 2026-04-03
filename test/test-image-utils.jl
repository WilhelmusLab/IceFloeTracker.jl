@testitem "apply_mask" begin
    using Images: RGB, Gray, N0f8

    @testset "apply_mask - Float64 matrix" begin
        img = [1.0 2.0; 3.0 4.0]
        mask = BitMatrix([true false; false true])
        result = apply_mask(img, mask)
        @test result == [0.0 2.0; 3.0 0.0]
        # original is unchanged
        @test img == [1.0 2.0; 3.0 4.0]
    end

    @testset "apply_mask - RGB image" begin
        img = [RGB(1.0, 0.5, 0.0) RGB(0.5, 1.0, 0.0); RGB(0.0, 0.5, 1.0) RGB(0.5, 0.0, 1.0)]
        mask = BitMatrix([true false; false true])
        result = apply_mask(img, mask)
        @test result[1, 1] == RGB(0.0, 0.0, 0.0)
        @test result[1, 2] == RGB(0.5, 1.0, 0.0)
        @test result[2, 1] == RGB(0.0, 0.5, 1.0)
        @test result[2, 2] == RGB(0.0, 0.0, 0.0)
        # original is unchanged
        @test img[1, 1] == RGB(1.0, 0.5, 0.0)
    end

    @testset "apply_mask - Gray image" begin
        img = Gray{Float64}.([0.5 0.8; 0.3 0.6])
        mask = BitMatrix([true false; false true])
        result = apply_mask(img, mask)
        @test result[1, 1] == Gray(0.0)
        @test result[1, 2] == Gray(0.8)
        @test result[2, 1] == Gray(0.3)
        @test result[2, 2] == Gray(0.0)
    end

    @testset "apply_mask - BitMatrix image" begin
        img = BitMatrix([true true; false true])
        mask = BitMatrix([true false; false true])
        result = apply_mask(img, mask)
        @test result == BitMatrix([false true; false false])
    end

    @testset "apply_mask! - in-place" begin
        img = [1.0 2.0; 3.0 4.0]
        mask = BitMatrix([true false; false true])
        apply_mask!(img, mask)
        @test img == [0.0 2.0; 3.0 0.0]
    end

    @testset "apply_mask! - allocates less than apply_mask" begin
        # Wrap in functions so compilation happens before measurement
        _test_copy(img, mask) = apply_mask(img, mask)
        _test_inplace!(img, mask) = apply_mask!(img, mask)
        img_copy = rand(100, 100)
        img_inplace = rand(100, 100)
        mask = BitMatrix(rand(Bool, 100, 100))
        # Warm up to avoid measuring compilation
        _test_copy(img_copy, mask)
        _test_inplace!(img_inplace, mask)
        x = @allocated _test_copy(img_copy, mask)
        y = @allocated _test_inplace!(img_inplace, mask)
        @test x > y
    end

    @testset "apply_mask - all false mask leaves image unchanged" begin
        img = [1.0 2.0; 3.0 4.0]
        mask = BitMatrix([false false; false false])
        result = apply_mask(img, mask)
        @test result == img
    end

    @testset "apply_mask - all true mask zeros image" begin
        img = [1.0 2.0; 3.0 4.0]
        mask = BitMatrix([true true; true true])
        result = apply_mask(img, mask)
        @test all(iszero, result)
    end
end
