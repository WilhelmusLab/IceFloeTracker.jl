@testset "overcast" begin
    imgdir = joinpath(test_data_dir, "pipeline/overcast")
    imgs = deserialize(joinpath(imgdir, "overcast.jls"))

    serialize(joinpath(imgdir, "generated_landmask.jls"), imgs.l)

    truecolordir = mkpath(joinpath(imgdir, "truecolor"))
    refdir = mkpath(joinpath(imgdir, "reflectance"))
    save(joinpath(truecolordir, "truecolor.png"), imgs.t)
    save(joinpath(refdir, "reflectance.png"), imgs.r)

    IceFloeTracker.Pipeline.preprocess(;
        truedir=truecolordir, refdir=refdir, lmdir=imgdir, output=imgdir
    )

    segmented_floes = deserialize(joinpath(imgdir, "segmented_floes.jls"))
    @test 0 == IceFloeTracker.nrow(
        IceFloeTracker.Pipeline.extractfeatures(segmented_floes[1]; features=["area"])
    )
    @test sum(segmented_floes[1]) == 0

    # clean up and restore original data
    rm(imgdir; recursive=true)
    mkpath(imgdir)
    serialize(joinpath(imgdir, "overcast.jls"), imgs)
end
