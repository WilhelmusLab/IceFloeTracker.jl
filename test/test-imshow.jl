@testitem "imshow Bitmatrix" begin
    import IceFloeTracker: imshow
    import Images: Gray

    bitmatrix = BitMatrix(rand(Bool, 10, 10))
    @test imshow(bitmatrix) isa Array{Gray{Bool},2}
end

@testitem "imshow fixed-point images" begin
    import IceFloeTracker: imshow
    import Images: N0f8, N0f16, N0f32, N0f64, Float32, Float64, Gray

    for T in (N0f8, N0f16, N0f32, N0f64, Float32, Float64)
        intmatrix = rand(T, 10, 10)
        @test imshow(intmatrix) isa Array{Gray{T},2}
    end
end
