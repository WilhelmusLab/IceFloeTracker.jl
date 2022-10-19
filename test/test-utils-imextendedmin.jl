@testset "imextendedmin" begin
    println("-------------------------------------------------")
    println("---------------- imextendedmin Tests ------------------")

    test_matrix = "$(test_data_dir)/test_extendedmin.csv"
    test_image = DelimitedFiles.readdlm(test_matrix, ',', Bool)
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
    matlab_extendedmin_output = DelimitedFiles.readdlm(
        matlab_extendedmin_output_file, ',', Bool
    )

    extendedmin_bitmatrix = IceFloeTracker.imextendedmin(test_image)
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
