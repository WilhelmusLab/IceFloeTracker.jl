@testitem "Watkins2026" begin
    using IceFloeTracker.Watkins2026:
        Dataset,
        Case,
        metadata,
        cases,
        modis_truecolor,
        modis_falsecolor,
        modis_landmask,
        modis_cloudfraction,
        masie_landmask,
        masie_seaice,
        validated_binary_floes,
        validated_labeled_floes,
        validated_floe_properties

    using Images: RGBA, N0f8, Colorant, SegmentedImage, Gray
    using DataFrames: nrow, DataFrame, DataFrameRow, subset

    dataset = Dataset()

    @testset "Watkins2026 Dataset" begin
        # Dataset interface
        @test dataset isa Dataset
        @test length(dataset) == 378
        @test nrow(metadata(dataset)) == 378

        # Cases interface
        cases_ = cases(dataset)
        @test length(cases_) == 378
        @info(cases_[1])
        @test cases_[1] isa Case
        @test cases_[end] isa Case
        @test metadata(cases_[1]) isa DataFrameRow
        @test metadata(cases_) isa DataFrame
    end

    @testset "Filtering and Subsetting" begin
        # Filtering interface
        filtered_dataset = filter(c -> c.case_number in (1, 2), dataset)
        @test filtered_dataset isa Dataset
        @test length(filtered_dataset) == 4
        @info filtered_dataset

        subsetted_dataset = subset(dataset, :case_number => c -> c .<= 2)
        @test subsetted_dataset isa Dataset
        @test length(subsetted_dataset) == 4
        @info subsetted_dataset
    end

    @testset "Watkins2026 Case Data" begin
        case = first(dataset)
        @test modis_truecolor(case) isa AbstractArray{RGBA{N0f8},2}
        @test modis_falsecolor(case) isa AbstractArray{RGBA{N0f8},2}
        @test modis_landmask(case) isa AbstractArray{<:Gray{Bool},2}
        @test modis_cloudfraction(case) isa AbstractArray{RGBA{N0f8},2}
        @test validated_binary_floes(case) isa AbstractArray{<:Gray{Bool},2}
        @test validated_labeled_floes(case) isa SegmentedImage
        @test validated_floe_properties(case) isa DataFrame
        @test masie_seaice(case) isa AbstractArray{<:Gray,2}
        @test masie_landmask(case) isa AbstractArray{<:Gray,2}
    end
end
