@testset "bwtraceboundary test" begin
    println("-------------------------------------------------")
    println("------------ bwtraceboundary Tests --------------")

    # Create an image with 3 connected components. The test consists of identifying the three closed sequences of border pixels in the image below. We do so using bwtraceboundary.
    A = zeros(Int, 13, 16)
    A[2:6, 2:6] .= 1
    A[4:8, 7:10] .= 1
    A[10:12, 13:15] .= 1
    A[10:12, 3:6] .= 1
    # 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    # 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0
    # 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0
    # 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
    # 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
    # 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
    # 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0
    # 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0
    # 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    # 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
    # 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
    # 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
    # 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

    # get boundaries closed by default, no point provided
    boundary = IceFloeTracker.bwtraceboundary(A)

    # Test 1: Check correct number of boundary pixels are obtained
    @test all([
        length(boundary[1]) == 27, length(boundary[2]) == 11, length(boundary[3]) == 9
    ])

    # Test 2: Initial points for contours
    p1 = (2, 2)
    p2 = (12, 3)
    p3 = (10, 15)
    pbad = (0, 0)
    pint = (4, 4)

    # get a closed boundary starting at p1
    out = IceFloeTracker.bwtraceboundary(A; P0=p1)
    @test all([length(boundary[1]) == length(out), out[1] == out[end]])

    # get a closed boundary starting at p2
    out = IceFloeTracker.bwtraceboundary(A; P0=p2)
    @test all([length(boundary[2]) == length(out), out[1] == out[end]])

    # get a closed boundary starting at p3
    out = IceFloeTracker.bwtraceboundary(A; P0=p3)
    @test all([length(boundary[3]) == length(out), out[1] == out[end]])

    # test exterior point
    out = IceFloeTracker.bwtraceboundary(A; P0=pbad)
    @test boundary == out

    # test interior point
    out = IceFloeTracker.bwtraceboundary(A; P0=pint)
    @test boundary == out

    # test input is a BitMatrix
    @test IceFloeTracker.bwtraceboundary(A) == IceFloeTracker.bwtraceboundary(BitArray(A))

    # test not closed
    @test length(IceFloeTracker.bwtraceboundary(A; P0=p1, closed=false)) ==
        length(IceFloeTracker.bwtraceboundary(A; P0=p1)) - 1
end;
