@testitem "hbreak tests" begin
    h1, h2 = keys(IceFloeTracker.make_hbreak_dict())
    A = zeros(Bool, 20, 20)

    # Make 4 H-connected groups on each corner of A
    A[1:3, 1:3] = h1
    A[1:3, (end - 2):end] = h2
    A[(end - 2):end, 1:3] = h1
    A[(end - 2):end, (end - 2):end] = h2
    A[10:12, 10:12] = h1 # and another in the middle for good measure

    @test sum(A) == sum(IceFloeTracker.hbreak(A)) + 5
end
