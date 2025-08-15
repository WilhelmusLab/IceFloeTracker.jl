@ntestset "find_ice_Labels" begin
    @ntestset "interface checks" begin
        @ntestset "functor version" begin
            data_loader = Watkins2025GitHub(;
                ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70"
            )
            case = first(data_loader(c -> (c.case_number == 12 && c.satellite == "terra")))
            landmask = case.modis_landmask
            falsecolor = case.modis_falsecolor
            algorithm = LopezAcosta2019IceDetection()
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
                    masker(.!(landmask))(falsecolor_image), LopezAcosta2019IceDetection()
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
                        LopezAcosta2019IceDetection(),
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
            find_ice(falsecolor, LopezAcosta2019IceDetection())
            find_ice(masker(landmask)(falsecolor), LopezAcosta2019IceDetection())
        end
    end
end
