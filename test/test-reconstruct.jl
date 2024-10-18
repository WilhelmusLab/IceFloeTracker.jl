# WIP
using ZipFile # delete later
using DelimitedFiles: readdlm
using IceFloeTracker: imcomplement

@testset "imcomplement" begin
    r = ZipFile.Reader("test_inputs/coins.zip")
    coins = readdlm(r.files[1], ',', Int)
    close(r)
    @test true


end
