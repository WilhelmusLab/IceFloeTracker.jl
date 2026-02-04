@testitem "k-means-segmentation" begin
    using IceFloeTracker
    import Images: Gray, segment_labels
    dataset = Watkins2026Dataset(; ref="v0.1")
    case = first(filter(c -> (c.case_number == 6 && c.satellite == "aqua"), dataset))
    
    ### vanilla case
    img = Gray.(modis_truecolor(case))
    kseg = kmeans_segmentation(img)
    @test length(segment_labels(kseg)) == 4

    kseg = kmeans_segmentation(img; k=5)
    @test length(segment_labels(kseg)) == 5

    ### tiled k-means
    tiles = get_tiles(img, 200)
    kseg_tiles = kmeans_segmentation(img, tiles)
    @test length(segment_labels(kseg_tiles) > 4)
end