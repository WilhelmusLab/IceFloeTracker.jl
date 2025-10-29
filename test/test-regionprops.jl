@testitem "regionprops" begin
    using Random
    import IceFloeTracker.Segmentation: regionprops
    import DataFrames: DataFrame, nrow
    import Images: label_components

    Random.seed!(123)
    bw_img = Bool.(rand([0, 1], 5, 10))
    bw_img[end, 7] = 1
    label_img = IceFloeTracker.label_components(bw_img, trues(3, 3))
    properties = (
        "centroid",
        "area",
        "major_axis_length",
        "minor_axis_length",
        "convex_area",
        "bbox",
        "perimeter",
        "orientation",
    )
    extra_props = nothing
    table = IceFloeTracker.regionprops_table(label_img, bw_img; properties=properties)
    total_labels = maximum(label_img)

    # Tests for regionprops_table

    @test typeof(table) <: DataFrame && # check correct data type
        6 == sum([p in names(table) for p in properties]) && # check correct set of properties
        size(table) == (total_labels, length(properties) + 4) # check correct table size

    # Check no value in bbox cols is 0 (zero indexing from skimage)
    @test all(Matrix(table[:, ["min_row", "min_col", "max_row", "max_col"]] .> 0))

    # Check no value in centroid cols is 0 (zero indexing from skimage)
    @test all(Matrix(table[:, ["row_centroid", "col_centroid"]] .> 0))

    # check default properties
    @test table == IceFloeTracker.regionprops_table(label_img)

    # Tests for regionprops
    regions = regionprops(label_img, bw_img)

    # Check some data matches in table
    randnum = rand(1:total_labels)
    @test table.area[randnum] == regions[randnum].area

    # Check floe masks generation and correct cropping
    IceFloeTracker.addfloemasks!(table, bw_img)
    @test all(
        [
            length(unique(label_components(table.mask[i], trues(3, 3)))) for
            i in 1:nrow(table)
        ] .== [2, 2, 1],
    )
end
