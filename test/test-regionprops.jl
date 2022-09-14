@testset "regionprops.jl" begin
    println("-------------------------------------------------")
    println("-------------- regionprops Tests ----------------")

    Random.seed!(123);
    bw_img = rand([0, 1], 5, 10);
    label_img = Images.label_components(bw_img, trues(3,3));
    properties = ["area", "perimeter", "area_convex","solidity"];
    extra_props = nothing;
    table = IceFloeTracker.regionprops_table(label_img, bw_img, properties = properties, dataframe=true)
    total_labels = maximum(label_img)
    
    # Tests for regionprops_table

    @test typeof(table) <: DataFrame && # check correct data type
          sort(names(table)) == sort(properties) && # check correct set of properties
          size(table) == (total_labels, length(properties)) # check correct table size

    dict = IceFloeTracker.regionprops_table(label_img, bw_img, properties = properties)
    @test typeof(dict) <: Dict && # check correct data type
          sort(collect(keys(dict))) == sort(properties) #  check correct size and properties

    # Tests for regionprops

    regions = IceFloeTracker.regionprops(label_img, bw_img);

    # Check some data matches in table
    randnum = rand(1:total_labels)
    @test table.area[randnum] == regions[randnum].area

end;
