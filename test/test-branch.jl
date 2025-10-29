@testitem "branch points tests" begin
    include("config.jl")
    import DelimitedFiles: readdlm
    import Images: erode

    dir = joinpath(test_data_dir, "branch")
    readcsv(f) = readdlm(joinpath(dir, f), ',', Bool)

    # Get inputs
    circles = readcsv("circles.csv")
    circles_skel = readcsv("circles_skel.csv")
    circles_branch_exp = readcsv("circles_branch_matlab.csv")

    # Ideal test on skeletonized image
    @test circles_branch_exp == branch(circles_skel)

    # Test on non-skeletonized image: Effect of brach = eroding
    @test erode(circles) == branch(circles)
end
