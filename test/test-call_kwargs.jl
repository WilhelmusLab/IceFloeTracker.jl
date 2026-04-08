@testitem "call_kwargs basic callback" begin
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

@testitem "call_kwargs with multiple functions" begin
    using IceFloeTracker.Utils: call_kwargs

    mutable1 = []
    mutable2 = []
    function func1(x)
        push!(mutable1, x)
        return nothing
    end
    function func2(x)
        push!(mutable2, x * 2)
        return nothing
    end

    callback = call_kwargs(; test1=func1, test2=func2)

    result = callback(; test1=3, test2=4)

    @test mutable1 == [3]
    @test mutable2 == [8]
end

@testitem "call_kwargs with println side-effect" begin
    using IceFloeTracker.Utils: call_kwargs

    output = IOBuffer()
    function save_to_disk(x)
        println(output, "Saving to disk: ", x)
        return nothing
    end

    callback = call_kwargs(; segmentation_result=save_to_disk)

    callback(; segmentation_result="[Segmented image data]")

    @test String(take!(output)) == "Saving to disk: [Segmented image data]\n"
end

@testitem "call_kwargs with save to disk side-effect" begin
    using IceFloeTracker.Utils: call_kwargs
    dir = tempdir()

    function save_to_disk(x)
        file = joinpath(dir, "result.txt")
        open(file, "w") do io
            println(io, "Saving to disk: ", x)
        end
        return nothing
    end

    callback = call_kwargs(; message=save_to_disk)

    result = callback(; message="[Data to write]")

    # Check that the file was created and contains the expected content
    file = joinpath(dir, "result.txt")  # This should be the same file name generated in save_to_disk
    @test isfile(file)

    content = read(file, String)
    @test content == "Saving to disk: [Data to write]\n"
end