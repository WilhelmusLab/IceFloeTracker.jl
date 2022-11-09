@testset "bridge tests" begin
    println("-------------------------------------------------")
    println("---------------- bridge tests -------------------")
    bwin = readdlm("./test_inputs/bridge/bridge_in.csv", ',', Bool)
    bwexpected = readdlm("./test_inputs/bridge/bridge_expected.csv", ',', Bool)
    @test bridge(bwin) == bwexpected
end
