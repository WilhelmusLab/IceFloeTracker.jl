@testset "persist.jl" begin
    println("-------------------------------------------------")
    println("---------- Persist Image Tests ------------")
    img_path = "/landmask.tiff"
    outimage_path = "outimage1.tiff"
    img = Images.load.(test_data_dir * img_path)

    # Test filename in variable
    IceFloeTracker.@persist img outimage_path
    @test isfile(outimage_path)

    # Test filename as string literal
    IceFloeTracker.@persist img "outimage2.tiff"
    @test isfile("outimage2.tiff")

    # Test no-filename call. Default filename startswith 'persisted_mask-' 
    # First clear all files that start with this prefix, if any
    [rm(f) for f in readdir() if startswith(f, "persisted_mask-")]
    @assert length([f for f in readdir() if startswith(f, "persisted_mask-")]) == 0

    IceFloeTracker.@persist img
    @test length([f for f in readdir() if startswith(f, "persisted_mask-")]) == 1

    # clean up - part 1!
    rm(outimage_path)
    rm("outimage2.tiff")

    # Part 2
    IceFloeTracker.@persist identity(img) outimage_path
    IceFloeTracker.@persist identity(img) "outimage2.tiff"
    IceFloeTracker.@persist identity(img) # no file name given
    @test isfile(outimage_path)
    @test isfile("outimage2.tiff")
    @test length([f for f in readdir() if startswith(f, "persisted_mask-")]) == 2

    # Part 3
    # Test filename as Expr, such as "$(dir)$(fname).$(ext)"
    dir = "./"
    fname = "persistedimg"
    ext = "png"
    IceFloeTracker.@persist img "$(dir)$(fname).$(ext)"
    @test isfile("$(dir)$(fname).$(ext)")

    # Clean up
    rm(outimage_path)
    rm("outimage2.tiff")
    rm("$(dir)$(fname).$(ext)")
    [rm(f) for f in readdir() if startswith(f, "persisted_mask-")]
end
