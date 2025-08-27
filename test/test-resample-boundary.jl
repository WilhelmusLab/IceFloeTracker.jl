@testitem "resample_boundary test" begin
    println("-------------------------------------------------")
    println("------------ resample_boundary Tests --------------")

    # Create an image with 3 connected components. The test consists of identifying the three closed sequences of border pixels in the image below. We do so using bwtraceboundary.
    A = zeros(Int, 9, 11)
    A[2:6, 2:6] .= 1
    A[4:8, 7:10] .= 1
    # 0  0  0  0  0  0  0  0  0  0  0
    # 0  1  1  1  1  1  0  0  0  0  0
    # 0  1  1  1  1  1  0  0  0  0  0
    # 0  1  1  1  1  1  1  1  1  1  0
    # 0  1  1  1  1  1  1  1  1  1  0
    # 0  1  1  1  1  1  1  1  1  1  0
    # 0  0  0  0  0  0  1  1  1  1  0
    # 0  0  0  0  0  0  1  1  1  1  0
    # 0  0  0  0  0  0  0  0  0  0  0

    # get boundary of biggest blob in image
    boundary = IceFloeTracker.bwtraceboundary(A; P0=(2, 2))

    # get resampled set of boundary points
    resampled_boundary = IceFloeTracker.resample_boundary(boundary)

    # Test 0: Check data type of resampled_boundary is Matrix{Float64}
    @test typeof(resampled_boundary) <: Matrix{Float64}

    # Test 1: Check correct number of resampled pixels where obtained
    @test size(resampled_boundary)[1] == length(boundary) ÷ 2

    # Test 2: Check perimeters are comparable (less than 10% after regularization by half the points)
    # make dist function
    norm(u) = (u[1]^2 + u[2]^2)^0.5

    per1 = sum([norm(u) for u in (boundary[2:end] - boundary[1:(end - 1)])])
    difs2 = [
        norm(u) for
        u in eachrow(resampled_boundary[2:end, :] - resampled_boundary[1:(end - 1), :])
    ]
    per2 = sum(difs2)
    d = abs(per1 - per2) / per1
    @test d < 0.1

    # Test 3: Check distances between a pair of adjacent points is about the same (small standard deviation)
    IceFloeTracker.std(difs2) < 1.0
end;
