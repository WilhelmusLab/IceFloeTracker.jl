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

    using Images: RGBA, N0f8, Colorant
    using DataFrames: nrow, DataFrame, subset

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

    filtered_dataset = filter(c -> c.case_number in (1, 2), dataset)
    @test filtered_dataset isa Dataset
    @test length(filtered_dataset) == 4
    @info filtered_dataset

    subsetted_dataset = subset(dataset, :case_number => c -> c .<= 2)
    @test subsetted_dataset isa Dataset
    @test length(subsetted_dataset) == 4
    @info subsetted_dataset

    for image_function in [
        modis_truecolor,
        modis_falsecolor,
        modis_landmask,
        modis_cloudfraction,
        masie_landmask,
        masie_seaice,
        validated_binary_floes,
        validated_labeled_floes,
    ]
        case = first(dataset)
        image = image_function(case)
        @test image isa AbstractArray{<:Colorant}
    end

    for table_function in [validated_floe_properties]
        case = first(dataset)
        table = table_function(case)
        @test table isa DataFrame
    end

    for image_function in [
        modis_truecolor,
        modis_falsecolor,
        modis_landmask,
        modis_cloudfraction,
        masie_landmask,
        masie_seaice,
        validated_binary_floes,
        validated_labeled_floes,
    ]
        smaller_filtered_dataset = filter(c -> c.case_number in (1,), dataset)
        images = image_function.(smaller_filtered_dataset)
        @test length(images) == 2
        @test images[1] isa AbstractArray{<:Colorant}
        @test images[end] isa AbstractArray{<:Colorant}
    end

    truecolor_images = modis_truecolor.(filtered_dataset)
    @test length(truecolor_images) == 4
    @test truecolor_images[1] isa AbstractArray{RGBA{N0f8},2}
    @test truecolor_images[end] isa AbstractArray{RGBA{N0f8},2}

    # filtered_dataset = Dataset(c -> c.case_number in (1, 2))
    # @test filtered_dataset isa Dataset
    # @test length(filtered_dataset) == 4

    # for case in dataset
    # @info case
    #     # @test case isa Case
    #     # @test modis_truecolor(case) isa AbstractArray{RGBA{N0f8},2}
    #     # @test modis_falsecolor(case) isa AbstractArray{RGBA{N0f8},2}
    #     # @test metadata(case) isa AbstractDict
    # end

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
