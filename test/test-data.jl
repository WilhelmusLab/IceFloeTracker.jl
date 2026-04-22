@testitem "Data" begin
    using Images: RGBA, N0f8, SegmentedImage, Gray
    using DataFrames: nrow, DataFrame, DataFrameRow, subset
    using Dates: DateTime

    @testset "Watkins2026Dataset" begin
        dataset = Watkins2026Dataset(; ref="v0.2")
        @testset "Dataset Properties" begin
            @test dataset isa Dataset
            @test length(dataset) == 378
            @test nrow(info(dataset)) == 378
            @test info(dataset) isa DataFrame
        end
        @testset "Case" begin
            dataset = Watkins2026Dataset(; ref="v0.2")
            @test dataset[1] isa Case
            @test info(dataset[1]) isa DataFrameRow
        end
        @testset "Filtering and Subsetting" begin
            dataset = Watkins2026Dataset(; ref="v0.2")
            filtered_dataset = filter(c -> c.case_number in (1, 2), dataset)
            @test filtered_dataset isa Dataset
            @test length(filtered_dataset) == 4
            subsetted_dataset = subset(dataset, :case_number => c -> c .<= 2)
            @test subsetted_dataset isa Dataset
            @test length(subsetted_dataset) == 4
        end
        @testset "Case Data" begin
            dataset = Watkins2026Dataset(; ref="v0.2")
            case = first(dataset)
            @test name(case) isa String
            @test pass_time(case) isa DateTime
            @test modis_truecolor(case) isa AbstractArray{RGBA{N0f8},2}
            @test modis_falsecolor(case) isa AbstractArray{RGBA{N0f8},2}
            @test modis_landmask(case) isa AbstractArray{<:Gray{Bool},2}
            @test modis_cloudfraction(case) isa AbstractArray{RGBA{N0f8},2}
            @test validated_binary_floes(case) isa AbstractArray{<:Gray{Bool},2}
            @test validated_binary_landfast(case) isa AbstractArray{<:Gray{Bool},2}
            @test validated_labeled_floes(case) isa SegmentedImage
            @test validated_floe_properties(case) isa DataFrame
            @test masie_seaice(case) isa AbstractArray{<:Gray,2}
            @test masie_landmask(case) isa AbstractArray{<:Gray,2}
        end
    end
end

@testitem "Data loader retries and validates files" begin
    using IceFloeTracker.Data: _get_file
    using Downloads: RequestError, Response

    @testset "retries failed download and succeeds" begin
        temp = mktempdir()
        target = joinpath(temp, "nested", "file.bin")
        attempts = Ref(0)

        download_fn =
            (url, path) -> begin
                attempts[] += 1
                attempts[] < 2 && throw(
                    RequestError(
                        url,
                        -1,
                        "transient",
                        Response("", "", 500, "Internal Server Error", []),
                    ),
                )
                open(path, "w") do io
                    write(io, UInt8(0x01))
                end
            end

        file = _get_file(
            "https://example.invalid/file.bin",
            target;
            max_attempts=3,
            download_fn=download_fn,
        )
        @test file == target
        @test attempts[] == 2
        @test isfile(target)
    end

    @testset "A generic permanent error just fails" begin
        temp = mktempdir()
        target = joinpath(temp, "nested", "file.bin")
        attempts = Ref(0)

        download_fn = (url, path) -> begin
            attempts[] += 1
            throw(ArgumentError("permanent"))
        end

        @test_throws ArgumentError _get_file(
            "https://example.invalid/file.bin", target; download_fn=download_fn
        )

        @test attempts[] == 1
    end

    @testset "An HTTP 429 Too Many Requests error stops retries" begin
        temp = mktempdir()
        target = joinpath(temp, "nested", "file.bin")
        attempts = Ref(0)

        download_fn =
            (url, path) -> begin
                attempts[] += 1
                throw(
                    RequestError(
                        url,
                        429,
                        "Too Many Requests",
                        Response("", "", 429, "Too Many Requests", []),
                    ),
                )
            end

        @test_throws RequestError _get_file(
            "https://example.invalid/file.bin", target; download_fn=download_fn
        )
        @test attempts[] == 1
    end

    @testset "invalid cached file triggers re-download" begin
        temp = mktempdir()
        target = joinpath(temp, "cached.csv")
        write(target, "")
        attempts = Ref(0)

        download_fn = (url, path) -> begin
            attempts[] += 1
            open(path, "w") do io
                write(io, "x")
            end
        end

        file = _get_file(
            "https://example.invalid/cached.csv",
            target;
            max_attempts=2,
            download_fn=download_fn,
        )
        @test file == target
        @test attempts[] == 1
        @test filesize(target) > 0
    end
end
