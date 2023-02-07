@testset "regionprops.jl" begin
    println("-------------------------------------------------")
    println("-------------- regionprops Tests ----------------")

    Random.seed!(123)
    bw_img = rand([0, 1], 5, 10)
    label_img = IceFloeTracker.label_components(bw_img, trues(3, 3))
    properties = ("centroid", "area", "major_axis_length", "minor_axis_length", "convex_area", "bbox")
    extra_props = nothing
    table = IceFloeTracker.regionprops_table(
        label_img, bw_img; properties=properties
    )
    total_labels = maximum(label_img)

    # Tests for regionprops_table

    @test typeof(table) <: DataFrame && # check correct data type
          4 == [p in names(table) for p in properties] |> sum && # check correct set of properties
          size(table) == (total_labels, length(properties) + 4) # check correct table size

    # Check no value in columns bbox-* of table is 0 (zero indexing from skimage)
    @test Matrix(table[:, ["bbox-0", "bbox-1", "bbox-2", "bbox-3"]] .> 0) |> all

    # check default properties
    @test table == IceFloeTracker.regionprops_table(label_img) 

    # Tests for regionprops

    regions = IceFloeTracker.regionprops(label_img, bw_img)

    # Check some data matches in table
    randnum = rand(1:total_labels)
    @test table.area[randnum] == regions[randnum].area
end
