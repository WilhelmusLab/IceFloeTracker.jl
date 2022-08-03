# remove these later
using Test
using Images
using DataFrames
bw_img = Images.load("test/test_inputs/matlab_landmask.png")
# @testset "regionprops.jl" begin
    println("-------------------------------------------------")
    println("-------------- regionprops Tests ----------------")

    
    properties = ["area",
                  "centroid",
                  "axis_major_length",
                  "axis_minor_length",
                  "orientation",
                  "area_convex",
                  "solidity"]
    regions_table = regionprops(bw_img, properties=properties, connectivity=2)
# end
