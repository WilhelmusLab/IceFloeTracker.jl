@testitem "regionprops" begin
    using Random
    import DataFrames: DataFrame, nrow
    import Images: label_components

    Random.seed!(123)
    bw_img = Bool.(rand([0, 1], 5, 10))
    bw_img[end, 7] = 1
    label_img = label_components(bw_img, trues(3, 3))
    properties = (
        "label",
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
    table = regionprops_table(label_img, bw_img; properties=properties)
    total_labels = maximum(label_img)

    # Tests for regionprops_table
    @test typeof(table) <: DataFrame # check correct data type
    nrow(table) == (total_labels) # check correct table size
    @info names(table)
    @test all(
        [
            "area",
            "min_row",
            "min_col",
            "max_row",
            "max_col",
            "row_centroid",
            "col_centroid",
            "convex_area",
            "label",
            "major_axis_length",
            "minor_axis_length",
            "orientation",
            "perimeter",
        ] ⊆ names(table),
    )
    # check all requested properties are present

    # Check no value in bbox cols is 0 (zero indexing from skimage)
    @test all(Matrix(table[:, ["min_row", "min_col", "max_row", "max_col"]] .> 0))

    # Check no value in centroid cols is 0 (zero indexing from skimage)
    @test all(Matrix(table[:, ["row_centroid", "col_centroid"]] .> 0))

    # check default properties
    @test table == regionprops_table(label_img)

    # Tests for regionprops
    regions = regionprops(label_img, bw_img)

    # Check some data matches in table
    randnum = rand(1:total_labels)
    @test table.area[randnum] == regions[randnum].area

    # Check floe masks generation and correct cropping
    add_floemasks!(table, label_img)
    @test all(
        [
            length(unique(label_components(table.mask[i], trues(3, 3)))) for
            i in 1:nrow(table)
        ] .== [2, 2, 1],
    )
end
