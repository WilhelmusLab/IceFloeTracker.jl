@testitem "Persist image" begin
    outimage_path = "outimage1.tiff"
    img = ones(3, 3)

    # Test filename in variable
    @persist img outimage_path
    @test isfile(outimage_path)

    # Test filename as string literal
    @persist img "outimage2.tiff"
    @test isfile("outimage2.tiff")

    # Test no-filename call. Default filename startswith 'persisted_img-' 
    # First clear all files that start with this prefix, if any
    prefix = "persisted_img-"
    [rm(f) for f in readdir() if startswith(f, prefix)]

    @persist img
    @test length([f for f in readdir() if startswith(f, prefix)]) == 1

    # clean up - part 1!
    rm(outimage_path)
    rm("outimage2.tiff")
    [rm(f) for f in readdir() if startswith(f, prefix)]

    # Part 2
    # Test call expressions
    @persist identity(img) outimage_path
    @persist identity(img) "outimage2.tiff"
    @persist identity(img) # no file name given
    @test isfile(outimage_path)
    @test isfile("outimage2.tiff")
    @test length([f for f in readdir() if startswith(f, prefix)]) == 1

    # Part 3
    # Test filename as Expr, such as "$(dir)$(fname).$(ext)"

    fname = "persistedimg"
    ext = "png"
    @persist img "$(fname).$(ext)"
    @test isfile("$(fname).$(ext)")

    # Part 4
    # Test timestamp

    [rm(f) for f in readdir() if startswith(f, "foo")] # clear "foo*" files

    # Persist twice the same img with different ts
    @persist img "foo.png" true
    foos = [f for f in readdir() if startswith(f, "foo")]
    @test length(foos) == 1
    @test length(foos[1]) == 25 # ts adds 19 chars 

    # Clean up
    [rm(f) for f in readdir() if endswith(f, "png") || endswith(f, "tiff")]
end
