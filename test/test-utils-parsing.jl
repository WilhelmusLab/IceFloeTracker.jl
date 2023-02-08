@testset verbose=true "utils.jl for arg parsing" begin

    @testset "utils.jl parse_2tuple" begin
        println("-------------------------------------------------")
        println("--------------- Parse 2 Tuple Tests -------------")
        s = raw"(1   ,   3)" 
        @test (1, 3) == IceFloeTracker.parse_2tuple(s)
        s = raw"(1,3)"
        @test (1, 3) == IceFloeTracker.parse_2tuple(s)
        s = raw"(1 3)"
        @test_throws ArgumentError IceFloeTracker.parse_2tuple(s)
        s = raw"(1, 3, 5)"
        @test_throws MethodError IceFloeTracker.parse_2tuple(s)
        s = "1 2"
        @test_throws ArgumentError IceFloeTracker.parse_2tuple(s)
        s= "1, 2)"
        @test_throws ArgumentError IceFloeTracker.parse_2tuple(s)
        s = "(1, 2"
        @test_throws ArgumentError IceFloeTracker.parse_2tuple(s)
    end

    @testset "utils.jl check_2_tuple" begin
        println("-------------------------------------------------")
        println("--------------- Check_2Tuple Tests -------------")
        @test_throws ArgumentError IceFloeTracker.check_2_tuple((1, 1))
        @test_throws ArgumentError IceFloeTracker.check_2_tuple((2, 1))
        @test IceFloeTracker.check_2_tuple((1, 2)) == nothing
        @test_throws MethodError IceFloeTracker.check_2_tuple((1, 2, 3))
        @test_throws MethodError IceFloeTracker.check_2_tuple((1,))
    end
end

