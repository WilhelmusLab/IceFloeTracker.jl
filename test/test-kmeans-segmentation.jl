@testitem "k-means-segmentation" begin
    using IceFloeTracker
    import Images: Gray, segment_labels
    dataset = Watkins2026Dataset(; ref="v0.1")
    case = first(filter(c -> (c.case_number == 6 && c.satellite == "aqua"), dataset))
    
    ### standard
    img = Gray.(modis_truecolor(case))
    kseg = kmeans_segmentation(img)
    @test length(segment_labels(kseg)) == 4

    kseg = kmeans_segmentation(img; k=5)
    @test length(segment_labels(kseg)) == 5

    ### tiled
    tiles = get_tiles(img, 200)
    kseg_tiles = kmeans_segmentation(img, tiles)
    @test length(segment_labels(kseg_tiles)) > 4
end

@testitem "k-means-binarization" begin
    using IceFloeTracker
    import Images: Gray, segment_labels
    dataset = Watkins2026Dataset(; ref="v0.1")
    case = first(filter(c -> (c.case_number == 6 && c.satellite == "aqua"), dataset))
    
    ### standard
    img = Gray.(modis_truecolor(case))
    fc_img = modis_falsecolor(case)
    algo =  IceDetectionThresholdMODIS721(band_7_max=0.2, band_2_min=0.3, band_1_min=0.3)
    kbin = kmeans_binarization(img, fc_img; cluster_selection_algorithm=algo)
    @test 0.53 < (sum(kbin) / prod(size(img))) < 0.54

    ### tiled
    tiles = get_tiles(img, 200)
    kbin = kmeans_binarization(img, fc_img, tiles; cluster_selection_algorithm=algo)
    @test 0.51 < (sum(kbin) / prod(size(img))) < 0.52
end

