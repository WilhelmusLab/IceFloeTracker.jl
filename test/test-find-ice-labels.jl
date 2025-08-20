using Images: ARGB

@ntestset "IceDetectionAlgorithm" begin
    @ntestset "IceDetectionThresholdMODIS721" begin
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

    @ntestset "IceDetectionBrightnessPeaksMODIS721" begin
        f = IceDetectionBrightnessPeaksMODIS721(5 / 255, 75 / 255)

        data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
        case = first(data_loader(c -> (c.case_number == 12 && c.satellite == "terra")))

        # There is only one pixel in this image marked as peak
        single_floe = case.modis_falsecolor[225:260, 110:160]
        @test sum(binarize(single_floe, f)) == 1

        just_water = case.modis_falsecolor[1:50, 350:400]
        @test sum(binarize(just_water, f)) == 0

        masked_land = masker(case.modis_landmask, case.modis_falsecolor)[1:50, 1:50]
        @test sum(binarize(masked_land, f)) == 0 broken = true
    end

    @ntestset "IceDetectionLopezAcosta2019" begin
        fc = IceDetectionLopezAcosta2019()
        f1 = fc.algorithms[1]
        f2 = fc.algorithms[2]
        f3 = fc.algorithms[3]

        data_loader = Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")
        case = first(data_loader(c -> (c.case_number == 12 && c.satellite == "terra")))

        bc = binarize(case.modis_falsecolor, fc)
        individual_results = [binarize(case.modis_falsecolor, f) for f in fc.algorithms]
        # Regardless of the input, the output should be identical to one of the algorithms
        @test any(bc == bi for bi in individual_results)
    end
end

@ntestset "find_ice_labels" begin
    @ntestset "interface checks" begin
        @ntestset "functor version" begin
            data_loader = Watkins2025GitHub(;
                ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70"
            )
            case = first(data_loader(c -> (c.case_number == 12 && c.satellite == "terra")))
            landmask = case.modis_landmask
            falsecolor = case.modis_falsecolor
            algorithm = IceDetectionLopezAcosta2019()
            @test find_ice(falsecolor, algorithm) == algorithm(falsecolor)
        end
    end
    @ntestset "find_ice_labels" begin
        @ntestset "matlab comparison" begin
            falsecolor_image = float64.(load(falsecolor_test_image_file)[test_region...])
            landmask = convert(BitMatrix, load(current_landmask_file))
            ice_labels_matlab = DelimitedFiles.readdlm(
                "$(test_data_dir)/ice_labels_matlab.csv", ','
            )
            ice_labels_matlab = vec(ice_labels_matlab)

            @ntestset "example 1" begin
                @time ice_labels_julia = IceFloeTracker.find_ice_labels(
                    falsecolor_image, landmask
                )
                DelimitedFiles.writedlm("ice_labels_julia.csv", ice_labels_julia, ',')
                @test ice_labels_matlab == ice_labels_julia
            end
            @ntestset "example 2" begin
                @time ice_labels_ice_floe_region = IceFloeTracker.find_ice_labels(
                    falsecolor_image[ice_floe_test_region...],
                    landmask[ice_floe_test_region...],
                )
                DelimitedFiles.writedlm(
                    "ice_labels_floe_region.csv", ice_labels_ice_floe_region, ','
                )
                @test ice_labels_ice_floe_region == [84787, 107015]
            end
        end
    end
    @ntestset "find_ice" begin
        @ntestset "matlab comparison" begin
            @ntestset "example 1" begin
                falsecolor_image =
                    float64.(load(falsecolor_test_image_file)[test_region...])
                landmask = convert(BitMatrix, load(current_landmask_file))
                ice_labels_matlab = DelimitedFiles.readdlm(
                    "$(test_data_dir)/ice_labels_matlab.csv", ','
                )
                ice_labels_matlab = vec(ice_labels_matlab)
                ice_binary_new = IceFloeTracker.find_ice(
                    masker(.!(landmask))(falsecolor_image), IceDetectionLopezAcosta2019()
                )
                ice_labels_julia_new = IceFloeTracker.get_ice_labels(ice_binary_new)
                @test ice_labels_julia_new == ice_labels_matlab
            end
            @ntestset "example 2" begin
                falsecolor_image =
                    float64.(load(falsecolor_test_image_file)[test_region...])
                landmask = convert(BitMatrix, load(current_landmask_file))
                ice_labels_ice_floe_region_new = IceFloeTracker.get_ice_labels(
                    IceFloeTracker.find_ice(
                        masker(.!(landmask))(falsecolor_image)[ice_floe_test_region...],
                        IceDetectionLopezAcosta2019(),
                    ),
                )
                @test ice_labels_ice_floe_region_new == [84787, 107015]
            end
        end
        @ntestset "validated data" begin
            data_loader = Watkins2025GitHub(;
                ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70"
            )
            case = first(data_loader(c -> (c.case_number == 12 && c.satellite == "terra")))
            landmask = case.modis_landmask
            falsecolor = case.modis_falsecolor
            baseline = find_ice(falsecolor, IceDetectionLopezAcosta2019())
            baseline_mask = find_ice(
                masker(landmask)(falsecolor), IceDetectionLopezAcosta2019()
            )
            fc_masked = masker(landmask)(falsecolor)

            @ntestset "IceDetectionLopezAcosta2019 type invariant" begin
                algorithm = IceDetectionLopezAcosta2019()
                @test find_ice(n0f8.(falsecolor), algorithm) == baseline
                @test find_ice(float64.(falsecolor), algorithm) == baseline
                @test find_ice(n4f12.(falsecolor), algorithm) == baseline broken = true
                @test find_ice(n0f8.(fc_masked), algorithm) == baseline_mask
                @test find_ice(float64.(fc_masked), algorithm) == baseline_mask
                @test find_ice(n4f12.(fc_masked), algorithm) == baseline_mask broken = true
            end

            @ntestset "IceDetectionThresholdMODIS721 type invariant" begin
                algorithm = IceDetectionThresholdMODIS721(;
                    band_7_max=(5 / 255), band_2_min=(230 / 255), band_1_min=(240 / 255)
                )
                @test find_ice(n0f8.(falsecolor), algorithm) == baseline
                @test find_ice(float64.(falsecolor), algorithm) == baseline
                @test find_ice(n4f12.(falsecolor), algorithm) == baseline broken = true
                @test find_ice(n0f8.(fc_masked), algorithm) == baseline_mask
                @test find_ice(float64.(fc_masked), algorithm) == baseline_mask
                @test find_ice(n4f12.(fc_masked), algorithm) == baseline_mask broken = true
            end
        end
    end
end
