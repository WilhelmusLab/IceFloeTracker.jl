@testset "prune tests" begin
    prune = IceFloeTracker.prune

    bw = Bool.([
        0 0 1 0 0
        0 0 1 0 0
        1 1 1 1 1
        0 0 1 0 0
        0 0 1 0 0
    ])
    @test Bool.([
        0 0 0 0 0
        0 0 1 0 0
        0 1 1 1 0
        0 0 1 0 0
        0 0 0 0 0
    ]) == prune(bw)

    sk = Bool.([
        0 0 1 0 0
        0 0 1 0 0
        1 1 1 0 0
        0 0 1 0 0
        0 0 1 0 0
    ])
    @test Bool.([
        0 0 0 0 0
        0 0 1 0 0
        0 1 1 0 0
        0 0 1 0 0
        0 0 0 0 0
    ]) == prune(sk)

    sk = Bool.([
        0 0 0 0 0
        0 0 0 0 0
        1 1 1 0 0
        0 0 1 0 0
        0 0 1 0 0
    ])
    @test Bool.([
        0 0 0 0 0
        0 0 0 0 0
        0 1 0 0 0
        0 0 1 0 0
        0 0 0 0 0
    ]) == prune(sk)

    sk = Bool.([
        0 0 0 0 0
        0 0 0 0 0
        1 1 1 1 1
        0 0 0 0 0
        0 0 0 0 0
    ])
    @test zeros(Bool, 5, 5) == prune(sk)
end
