@testitem "imshow Bitmatrix" begin
    import IceFloeTracker: imshow
    import Images: Gray

    bitmatrix = BitMatrix(rand(Bool, 10, 10))
    @test imshow(bitmatrix) isa Array{Gray{Bool},2}
end

@testitem "imshow fixed-point images" begin
    import IceFloeTracker: imshow
    import Images: FixedPoint, N0f8, N0f16, N0f32, N0f64, N4f12, N4f28, Gray

    for T in [N0f8, N0f16, N0f32, N0f64, N4f12, N4f28, Float16, Float32, Float64, BigFloat]
        array = rand(T, 10, 10)
        @test imshow(array) isa Array{Gray{T},2}
    end
end

@testitem "imshow floats outside [0, 1]" begin
    import IceFloeTracker: imshow

    for T in [Float16, Float32, Float64]
        array = rand(T, 10, 10) .* 10 .- 5 # values in range [-5, 5]
        @test imshow(array) isa Array{Gray{T},2}
    end
end
