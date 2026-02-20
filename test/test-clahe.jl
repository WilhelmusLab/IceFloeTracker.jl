@testitem "Initialization tests" begin
    @test ContrastLimitedAdaptiveHistogramEqualization() isa
        ContrastLimitedAdaptiveHistogramEqualization
end

@testitem "Parameter validation tests" begin
    using IceFloeTracker.Filtering.CLAHE: validate_parameters
    @test_throws ArgumentError validate_parameters(
        ContrastLimitedAdaptiveHistogramEqualization(; rblocks=0)
    )
    @test_throws ArgumentError validate_parameters(
        ContrastLimitedAdaptiveHistogramEqualization(; rblocks=-1)
    )
    @test_throws ArgumentError validate_parameters(
        ContrastLimitedAdaptiveHistogramEqualization(; cblocks=0)
    )
    @test_throws ArgumentError validate_parameters(
        ContrastLimitedAdaptiveHistogramEqualization(; cblocks=-1)
    )
end

@testitem "No crash for different image sizes" begin
    using Images, TestImages
    f = ContrastLimitedAdaptiveHistogramEqualization()
    image = "cameraman" # Grayscale
    for xsz in [63, 64, 65, 96, 128, 256, 300, 512],
        ysz in [63, 64, 65, 96, 128, 256, 300, 512]

        @info "Testing size: ($xsz, $ysz) for image: $image"
        img = imresize(testimage(image), (xsz, ysz))
        out = adjust_histogram(img, f)
        @test size(out) == size(img)
    end
end

@testitem "No crash for different image types" begin
    using Images, TestImages
    f = ContrastLimitedAdaptiveHistogramEqualization()
    image = "cameraman"
    img = testimage(image)

    @test adjust_histogram(convert.(Gray{N0f8}, img), f) isa Array{Gray{N0f8},2}
    @test adjust_histogram(convert.(Gray{N0f16}, img), f) isa Array{Gray{N0f16},2}
    @test adjust_histogram(convert.(Gray{Float64}, img), f) isa Array{Gray{Float64},2}
    @test adjust_histogram(convert.(Gray{Float32}, img), f) isa Array{Gray{Float32},2}

    @test adjust_histogram(convert.(RGB{N0f8}, img), f) isa Array{RGB{N0f8},2}
    @test adjust_histogram(convert.(RGB{N0f16}, img), f) isa Array{RGB{N0f16},2}
    @test adjust_histogram(convert.(RGB{Float64}, img), f) isa Array{RGB{Float64},2}
    @test adjust_histogram(convert.(RGB{Float32}, img), f) isa Array{RGB{Float32},2}

    @test adjust_histogram(convert.(RGBA{N0f8}, img), f) isa Array{RGBA{N0f8},2} broken =
        true
    @test adjust_histogram(convert.(RGBA{N0f16}, img), f) isa Array{RGBA{N0f16},2} broken =
        true
    @test adjust_histogram(convert.(RGBA{Float64}, img), f) isa Array{RGBA{Float64},2} broken =
        true
    @test adjust_histogram(convert.(RGBA{Float32}, img), f) isa Array{RGBA{Float32},2} broken =
        true
end
