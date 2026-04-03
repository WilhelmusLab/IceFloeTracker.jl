@testitem "call_kwargs callback" begin
    using IceFloeTracker.Utils: call_kwargs

    mutable = []
    function test_func(x)
        push!(mutable, x)
        return nothing
    end

    callback = call_kwargs(; test=test_func)

    result = callback(; test=5)

    @test mutable == [5]
end