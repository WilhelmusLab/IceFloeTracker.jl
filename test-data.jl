@testitem "Data" begin
    using Images: RGBA, N0f8, SegmentedImage, Gray
    using DataFrames: nrow, DataFrame, DataFrameRow, subset
    using IceFloeTracker: Case, Dataset, metadata, loader, Watkins2026Dataset

    @testset "Watkins2026Dataset" begin
        dataset = Watkins2026Dataset(; ref="b865acc62f223d6ff14a073a297d682c4c034e5d")
        @testset "Dataset Properties" begin
            @test dataset isa Dataset
            @test length(dataset) == 378
            @test nrow(metadata(dataset)) == 378
            @test metadata(dataset) isa DataFrame
        end
        @testset "Case" begin
            dataset = Watkins2026Dataset(; ref="b865acc62f223d6ff14a073a297d682c4c034e5d")
            @test dataset[1] isa Case
            @test metadata(dataset[1]) isa DataFrameRow
        end
        @testset "Filtering and Subsetting" begin
            dataset = Watkins2026Dataset(; ref="b865acc62f223d6ff14a073a297d682c4c034e5d")
            filtered_dataset = filter(c -> c.case_number in (1, 2), dataset)
            @test filtered_dataset isa Dataset
            @test length(filtered_dataset) == 4
            subsetted_dataset = subset(dataset, :case_number => c -> c .<= 2)
            @test subsetted_dataset isa Dataset
            @test length(subsetted_dataset) == 4
        end
        @testset "Case Data" begin
            dataset = Watkins2026Dataset(; ref="b865acc62f223d6ff14a073a297d682c4c034e5d")
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
end
