@testitem "Watkins2026" begin
    using IceFloeTracker.Watkins2026: Dataset, Case, metadata, cases
    using DataFrames: nrow, DataFrame

    # Ideal interface:
    dataset = Dataset()
    @test dataset isa Dataset
    @test length(dataset) == 378
    @test nrow(metadata(dataset)) == 378

    cases_ = cases(dataset)
    @test length(cases_) == 378
    @info(cases_[1])
    @test cases_[1] isa Case
    @test cases_[end] isa Case
    @test metadata(cases_[1]) isa AbstractDict
    @test metadata(cases_) isa DataFrame

    # filtered_dataset = Dataset(c -> c.case_number in (1, 2))
    # @test filtered_dataset isa Dataset
    # @test length(filtered_dataset) == 4

    for case in dataset
        @info case
        #     # @test case isa Case
        #     # @test modis_truecolor(case) isa AbstractArray{RGBA{N0f8},2}
        #     # @test modis_falsecolor(case) isa AbstractArray{RGBA{N0f8},2}
        #     # @test metadata(case) isa AbstractDict
    end

    # @test nrow(metadata(loader)) == 180

    # dataset = data_loader(c -> c.case_number == 1)

    # @test dataset isa Dataset
    # @test metadata(dataset.metadata isa DataFrame
    # @test nrow(dataset.metadata) == 2
    # @info dataset.metadata

    #@test length(cases) == 2

    # @test length(modis_truecolor.(cases)) == 2
    # @test length(modis_falsecolor.(cases)) == 2

    # @test modis_truecolor(first(cases)) isa Array{RGB{N0f8},2}
end
