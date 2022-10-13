@testset "bwareamaxfilt test" begin
    println("-------------------------------------------------")
    println("------------ bwareamaxfilt Tests --------------")

    # Create a bitmatrix with a big floe and two smaller floes -- 3 connected components in total.
    A = zeros(Bool, 12, 15); A[2:6, 2:6] .= 1; A[4:8, 7:10] .= 1; A[10:12, 13:15] .= 1; A[10:12, 3:6] .= 1;
    #  12×15 Matrix{Bool}:
    #  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    #  0  1  1  1  1  1  0  0  0  0  0  0  0  0  0
    #  0  1  1  1  1  1  0  0  0  0  0  0  0  0  0
    #  0  1  1  1  1  1  1  1  1  1  0  0  0  0  0
    #  0  1  1  1  1  1  1  1  1  1  0  0  0  0  0
    #  0  1  1  1  1  1  1  1  1  1  0  0  0  0  0
    #  0  0  0  0  0  0  1  1  1  1  0  0  0  0  0
    #  0  0  0  0  0  0  1  1  1  1  0  0  0  0  0
    #  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    #  0  0  1  1  1  1  0  0  0  0  0  0  1  1  1
    #  0  0  1  1  1  1  0  0  0  0  0  0  1  1  1
    #  0  0  1  1  1  1  0  0  0  0  0  0  1  1  1

    # Test 1: Check number of total blobs and their respective areas.
    # First get the labels
    labels = label_components(A);
    # 12×15 Matrix{Int64}:
    # 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    # 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0
    # 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0
    # 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0
    # 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0
    # 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0
    # 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0
    # 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0
    # 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    # 0  0  2  2  2  2  0  0  0  0  0  0  3  3  3
    # 0  0  2  2  2  2  0  0  0  0  0  0  3  3  3
    # 0  0  2  2  2  2  0  0  0  0  0  0  3  3  3

    # a) Test there are three blobs
    d = IceFloeTracker.get_areas(labels); @test length(d) == 3

    # b) Test largest blob is the one with label '1'
    @test IceFloeTracker.get_max_label(d) == 1

    # c) Test the distribution of the labels
    @test all([d[1]==45,d[2] == 12, d[3] == 9])

    # Test 2: Filter smaller blobs from label matrix
    @test sum(IceFloeTracker.filt_except_label(labels, IceFloeTracker.get_max_label(d)) .!=0) == 45

    # Test 3: Keep largest blob in input matrix a
    @test sum(IceFloeTracker.bwareamaxfilt(A)) == 45

    # Test 4: In-place version of bwareamaxfilt
    A_copy = copy(A);
    IceFloeTracker.bwareamaxfilt!(A_copy)
    @test A_copy == IceFloeTracker.bwareamaxfilt(A)
end
