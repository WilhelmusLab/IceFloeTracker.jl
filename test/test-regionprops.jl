@testitem "regionprops" begin
    using Random
    import DataFrames: DataFrame, nrow
    import Images: label_components

    Random.seed!(123)
    bw_img = Bool.(rand([0, 1], 5, 10))
    bw_img[end, 7] = 1
    label_img = label_components(bw_img, trues(3, 3))
    properties = [
        "label",
        "centroid",
        "area",
        "major_axis_length",
        "minor_axis_length",
        "bbox",
        "perimeter",
        "orientation"
    ]
    extra_props = nothing

    table = regionprops_table(label_img, bw_img; properties=properties, minimum_area=1)
    total_labels = maximum(label_img)
    
    # Tests for regionprops_table
    @test typeof(table) <: DataFrame # check correct data type
    @test nrow(table) == (total_labels - 1) # One point has area 1, we exclude that point intentionally
    @test all(
        [
            "area",
            "min_row",
            "min_col",
            "max_row",
            "max_col",
            "row_centroid",
            "col_centroid",
            "label",
            "major_axis_length",
            "minor_axis_length",
            "orientation",
            "perimeter",
        ] ⊆ names(table),
    )
    # check all requested properties are present

    # check default properties
    # TODO: verify that this test was needed. It doesn't seem useful to me.
    # @test table == regionprops_table(label_img)

    # Tests for regionprops
    regions = regionprops(label_img, bw_img)

    @test table.area[1] == regions[1].area

    # Check floe masks generation and correct cropping
    add_floemasks!(table, label_img)
    @test all(
        [
            length(unique(label_components(table.mask[i], trues(3, 3)))) for
            i in 1:nrow(table)
        ] .== [2, 2],
    )
end
