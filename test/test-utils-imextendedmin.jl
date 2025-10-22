@testitem "imextendedmin and bwdist" begin
    using DelimitedFiles: readdlm

    include("config.jl")

    test_matrix = "$(test_data_dir)/test_extendedmin.csv"
    test_image = readdlm(test_matrix, ',', Bool)
    # Test matrix
    # 10×10 BitMatrix:
    # 1 1 1 1 1 1 1 1 1 0 
    # 1 1 1 1 1 1 1 0 1 0
    # 1 1 1 0 1 1 1 1 1 1
    # 1 1 0 0 1 1 1 1 1 1
    # 1 1 1 1 1 1 1 1 1 1
    # 1 1 1 1 1 1 1 1 1 1
    # 1 1 1 1 1 0 1 1 1 1
    # 1 1 1 1 1 1 1 1 0 1
    # 1 0 1 1 1 1 1 1 1 1
    # 1 0 1 1 1 1 1 1 1 1

    matlab_extendedmin_output_file = "$(test_data_dir)/test_extendedmin_output.csv"
    matlab_extendedmin_output = readdlm(matlab_extendedmin_output_file, ',', Bool)

    # Test workflow for watershed segmentation
    distances = -IceFloeTracker.bwdist(.!test_image)
    extendedmin_bitmatrix = IceFloeTracker.imextendedmin(distances)
    # Matlab output
    # 10×10 BitMatrix:
    # 1 1 1 1 1 1 0 0 0 0
    # 1 1 0 0 0 1 0 0 0 0
    # 1 0 0 0 0 1 0 0 0 0
    # 1 0 0 0 0 1 1 1 1 1
    # 1 0 0 0 0 1 1 1 1 1
    # 1 1 1 1 0 0 0 1 1 1
    # 1 1 1 1 0 0 0 0 0 0
    # 0 0 0 1 0 0 0 0 0 0
    # 0 0 0 1 1 1 1 0 0 0
    # 0 0 0 1 1 1 1 1 1 1

    @test extendedmin_bitmatrix == matlab_extendedmin_output
end
