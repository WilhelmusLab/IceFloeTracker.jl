
@testitem "IceDetectionAlgorithm" begin
    using Images: ARGB, binarize, N0f8, RGB, Gray, load

    @testset "IceDetectionThresholdMODIS721" begin
        f = IceDetectionThresholdMODIS721(0.02, 0.92, 0.92)

        # Input colors
        water = ARGB{N0f8}(3 / 255, 13 / 255, 29 / 255, 1)
        sea_ice = ARGB{N0f8}(0 / 255, 237 / 255, 239 / 255, 1)
        grass = ARGB{N0f8}(96 / 255, 137 / 255, 71 / 255, 1)
        masked = ARGB{N0f8}(0 / 255, 237 / 255, 239 / 255, 0)

        shadow = RGB{N0f8}(0 / 255, 70 / 255, 117 / 255)
        land_ice = RGB{N0f8}(0 / 255, 237 / 255, 239 / 255)
        masked_land_ice = ARGB{N0f8}(0 / 255, 237 / 255, 239 / 255, 0)

        # Output colors
        is_ice = Gray(1)
        is_not_ice = Gray(0)

        @test binarize([water], f) == [is_not_ice]
        @test binarize([sea_ice], f) == [is_ice]
        @test binarize([grass], f) == [is_not_ice]
        @test binarize([masked], f) == [is_not_ice]
        @test binarize([shadow], f) == [is_not_ice]
        @test binarize([land_ice], f) == [is_ice]
        @test binarize([masked_land_ice], f) == [is_not_ice]
    end

    @testset "IceDetectionBrightnessPeaksMODIS721" begin
        f = IceDetectionBrightnessPeaksMODIS721(band_7_max = 5 / 255, possible_ice_threshold = 75 / 255)

        dataset = Watkins2026Dataset(; ref="v0.1")
        case = first(filter(c -> (c.case_number == 12 && c.satellite == "terra"), dataset))

        # There is only one pixel in this image marked as peak
        single_floe = modis_falsecolor(case)[225:260, 110:160]
        @test sum(binarize(single_floe, f)) == 77 # Changed from 1, new method picks up more ice

        just_water = modis_falsecolor(case)[1:50, 350:400]
        @test sum(binarize(just_water, f)) == 2 # Changed from 0, there are some ice pixels in the top of the image

        masked_land = masker(modis_landmask(case), modis_falsecolor(case))[1:50, 1:50]
        @test sum(binarize(masked_land, f)) == 0

        # Make sure we can all the new option
        f2 = IceDetectionBrightnessPeaksMODIS721(band_7_max = 5 / 255, possible_ice_threshold = 75 / 255, join_method="union")
        intersect_method = binarize(modis_falsecolor(case), f)
        union_method = binarize(modis_falsecolor(case), f2)
        @test sum(intersect_method) .<= sum(union_method)

        # Test whether it will default to intersect as intended
        f3 = IceDetectionBrightnessPeaksMODIS721(band_7_max = 5 / 255, possible_ice_threshold = 75 / 255, join_method="divide")
        alt_method = binarize(modis_falsecolor(case), f3)
        @test allequal(alt_method .== intersect_method)
    end

    @testset "IceDetectionLopezAcosta2019" begin
        fc = IceDetectionLopezAcosta2019()

        dataset = Watkins2026Dataset(; ref="v0.1")
        case = first(filter(c -> (c.case_number == 12 && c.satellite == "terra"), dataset))

        bc = binarize(modis_falsecolor(case), fc)
        individual_results = [binarize(modis_falsecolor(case), f) for f in fc.algorithms]
        # Regardless of the input, the output should be identical to one of the algorithms
        @test any(bc == bi for bi in individual_results)
        # The first non-zero output should be the same as the return value
        first_non_zero = first(bi for bi in individual_results if sum(bi) > 0)
        @test bc == first_non_zero
    end
end

@testitem "find_ice_labels" begin
    using Images: binarize, n0f8, float64, n4f12, n0f8, float64, n4f12, load, Gray
    import DelimitedFiles: readdlm, writedlm

    include("config.jl")

    @testset "interface checks" begin
        @testset "functor version" begin
            dataset = Watkins2026Dataset(; ref="v0.1")
            case = first(
                filter(c -> (c.case_number == 12 && c.satellite == "terra"), dataset)
            )
            landmask = modis_landmask(case)
            falsecolor = modis_falsecolor(case)
            algorithm = IceDetectionLopezAcosta2019()
            @test binarize(falsecolor, algorithm) == algorithm(falsecolor)
        end
    end
    @testset "find_ice_labels" begin
        @testset "matlab comparison" begin
            falsecolor_image = float64.(load(falsecolor_test_image_file)[test_region...])
            landmask = convert(BitMatrix, load(current_landmask_file)[test_region...])
            ice_labels_matlab = readdlm("$(test_data_dir)/ice_labels_matlab.csv", ',')
            ice_labels_matlab = vec(ice_labels_matlab)

            @testset "example 1" begin
                @time ice_labels_julia = find_ice_labels(falsecolor_image, landmask)
                writedlm("$(test_output_dir)/ice_labels_julia.csv", ice_labels_julia, ',')
                @test ice_labels_matlab == ice_labels_julia
            end
            @testset "example 2" begin
                @time ice_labels_ice_floe_region = find_ice_labels(
                    falsecolor_image[ice_floe_test_region...],
                    landmask[ice_floe_test_region...],
                )
                writedlm(
                    "$(test_output_dir)/ice_labels_floe_region.csv",
                    ice_labels_ice_floe_region,
                    ',',
                )
                @test ice_labels_ice_floe_region == [84787, 107015]
            end
        end
        @testset "get_ice_peaks" begin
            using Random
            using Images: build_histogram
            Random.seed!(123)
            img = Gray.(rand(0:255, 10, 10) ./ 255)
            edges, counts = build_histogram(img, 64; minval=0, maxval=1)
            pk = get_ice_peaks(edges, counts)
            @test pk == 0.375

            # Counts from case "111-greenland_sea-100km-20120623-terra-250m"
            edges = 0.0:0.015625:0.984375
            counts = [
                0.22849,
                0.04548,
                0.0317,
                0.02294,
                0.01893,
                0.01466,
                0.01293,
                0.011,
                0.00977,
                0.00934,
                0.00841,
                0.00838,
                0.00746,
                0.00743,
                0.00736,
                0.00682,
                0.00684,
                0.00667,
                0.00638,
                0.00712,
                0.00615,
                0.0065,
                0.00614,
                0.00601,
                0.00612,
                0.00683,
                0.00734,
                0.00782,
                0.00882,
                0.0092,
                0.00828,
                0.00926,
                0.00923,
                0.00933,
                0.00934,
                0.00978,
                0.0096,
                0.01093,
                0.01162,
                0.01295,
                0.0149,
                0.01803,
                0.02223,
                0.0239,
                0.02755,
                0.03081,
                0.03796,
                0.04513,
                0.054,
                0.06292,
                0.06747,
                0.07155,
                0.08288,
                0.07175,
                0.09395,
                0.0995,
                0.00321,
                0.00075,
                0.00022,
                4.0e-5,
                1.0e-5,
                0.0,
                0.0,
                0.0,
            ]
            pk = get_ice_peaks(edges, counts)
            @test pk == 0.859375
        end
    end
    @testset "binarize" begin
        @testset "matlab comparison" begin
            @testset "example 1" begin
                falsecolor_image =
                    float64.(load(falsecolor_test_image_file)[test_region...])
                landmask = convert(BitMatrix, load(current_landmask_file)[test_region...])
                ice_labels_matlab = readdlm("$(test_data_dir)/ice_labels_matlab.csv", ',')
                ice_labels_matlab = vec(ice_labels_matlab)
                ice_binary_new = binarize(
                    masker(landmask)(falsecolor_image), IceDetectionLopezAcosta2019()
                )
                ice_labels_julia_new = get_ice_labels(ice_binary_new)
                @test ice_labels_julia_new == ice_labels_matlab
            end
            @testset "example 2" begin
                falsecolor_image =
                    float64.(load(falsecolor_test_image_file)[test_region...])
                landmask = convert(BitMatrix, load(current_landmask_file)[test_region...])
                ice_labels_ice_floe_region_new = get_ice_labels(
                    binarize(
                        masker(landmask)(falsecolor_image)[ice_floe_test_region...],
                        IceDetectionLopezAcosta2019(),
                    ),
                )
                @test ice_labels_ice_floe_region_new == [84787, 107015]
            end
        end
        @testset "validated data" begin
            dataset = Watkins2026Dataset(; ref="v0.1")
            case = first(
                filter(c -> (c.case_number == 12 && c.satellite == "terra"), dataset)
            )
            landmask = modis_landmask(case)
            falsecolor = modis_falsecolor(case)
            baseline = binarize(falsecolor, IceDetectionLopezAcosta2019())
            baseline_mask = binarize(
                masker(landmask)(falsecolor), IceDetectionLopezAcosta2019()
            )
            fc_masked = masker(landmask)(falsecolor)

            @testset "IceDetectionLopezAcosta2019 type invariant" begin
                algorithm = IceDetectionLopezAcosta2019()
                @test binarize(n0f8.(falsecolor), algorithm) == baseline
                @test binarize(float64.(falsecolor), algorithm) == baseline
                @test binarize(n4f12.(falsecolor), algorithm) == baseline broken = true
                @test binarize(n0f8.(fc_masked), algorithm) == baseline_mask
                @test binarize(float64.(fc_masked), algorithm) == baseline_mask
                @test binarize(n4f12.(fc_masked), algorithm) == baseline_mask broken = true
            end

            @testset "IceDetectionThresholdMODIS721 type invariant" begin
                algorithm = IceDetectionThresholdMODIS721(;
                    band_7_max=(5 / 255), band_2_min=(230 / 255), band_1_min=(240 / 255)
                )
                @test binarize(n0f8.(falsecolor), algorithm) == baseline
                @test binarize(float64.(falsecolor), algorithm) == baseline
                @test binarize(n4f12.(falsecolor), algorithm) == baseline broken = true
                @test binarize(n0f8.(fc_masked), algorithm) == baseline_mask
                @test binarize(float64.(fc_masked), algorithm) == baseline_mask
                @test binarize(n4f12.(fc_masked), algorithm) == baseline_mask broken = true
            end
        end
    end
end

