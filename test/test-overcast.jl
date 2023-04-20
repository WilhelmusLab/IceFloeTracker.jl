@testset "overcast" begin
    imgdir = joinpath(test_data_dir, "pipeline/overcast")
    mkpath(imgdir)
    r = convert(Matrix{RGB{Float64}}, rand(10, 10))
    t = convert(Matrix{RGB{Float64}}, rand(10, 10))
    c = ones(Bool, 10, 10) # fully occluded
    l = (dilated=rand(Bool, 10, 10), non_dilated=ones(Bool, 10, 10))

    imgs = (t=truecolor_image, r=reflectance_image, l=landmask_imgs, c=cloudmask)

    serialize(joinpath(imgdir, "generated_landmask.jls"), imgs.l)

    truecolordir = mkpath(joinpath(imgdir, "truecolor"))
    refdir = mkpath(joinpath(imgdir, "reflectance"))
    save(joinpath(truecolordir, "truecolor.png"), imgs.t)
    save(joinpath(refdir, "reflectance.png"), imgs.r)

    IceFloeTracker.Pipeline.preprocess(;
        truedir=truecolordir, refdir=refdir, lmdir=imgdir, output=imgdir
    )

    segmented_floes = deserialize(joinpath(imgdir, "segmented_floes.jls"))

    # check output img and corresponding props df are empty
    @test 0 == IceFloeTracker.nrow(
        IceFloeTracker.Pipeline.extractfeatures(segmented_floes[1]; features=["area"])
    ) # empty props df
    @test sum(segmented_floes[1]) == 0 # empty image

    # clean up
    rm(imgdir; recursive=true)
end
