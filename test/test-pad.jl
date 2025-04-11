using IceFloeTracker: pad_images

@testset "padding" begin
    @testset "two identical inputs are padded the same" begin
        function test_padding_for_identical_inputs(image, expected_result)
            # display(image)
            result = pad_images(image, image)
            # display(result[1])
            @test result[1] == result[2]  # results are the same
            @test result[1] == expected_result  # results are as expected
        end

        @testset "1x1 zeros" begin
            test_padding_for_identical_inputs(
                Bool[0;;],
                Bool[
                    0 0 0
                    0 0 0
                    0 0 0
                ]
            )
        end

        @testset "3x3 one centered" begin
            test_padding_for_identical_inputs(
                Bool[1;;],
                Bool[
                    0 0 0
                    0 1 0
                    0 0 0
                ]
            )
        end

        @testset "3x3 ones" begin
            test_padding_for_identical_inputs(
                Bool[
                    1 1 1
                    1 1 1
                    1 1 1
                ],
                Bool[
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 1 1 1 0 0 0
                    0 0 0 1 1 1 0 0 0
                    0 0 0 1 1 1 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                ]
            )
        end

        @testset "2x2 anti-diagonal" begin
            test_padding_for_identical_inputs(
                Bool[
                    0 1
                    1 0
                ],
                Bool[
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                    0 0 0 1 0 0
                    0 0 1 0 0 0
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                ]
            )
        end
        @testset "2x2 diagonal" begin
            test_padding_for_identical_inputs(
                Bool[
                    1 0
                    0 1
                ],
                Bool[
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                    0 0 1 0 0 0
                    0 0 0 1 0 0
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                ]
            )
        end
        @testset "1x2 ones" begin
            test_padding_for_identical_inputs(
                Bool[
                    1 1;
                ],
                Bool[
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                    0 0 1 1 0 0
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                ]
            )
        end
        @testset "2x1 ones" begin
            test_padding_for_identical_inputs(
                Bool[
                    1;
                    1;;
                ],
                Bool[
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 1 0 0
                    0 0 1 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                ]
            )
        end
        @testset "1x3 ones" begin
            test_padding_for_identical_inputs(
                Bool[
                    1 1 1;
                ],
                Bool[
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 1 1 1 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                ]
            )
        end
    end
    @testset "unlike images are padded similarly" begin
        function test_padding_for_unlike_images(image1, image2, expected_result1, expected_result2)
            result1, result2 = pad_images(image1, image2)
            show(result1)
            @test result1 == expected_result1
            show(result2)
            @test result2 == expected_result2
        end
        @testset "minimum example" begin
            test_padding_for_unlike_images(
                Bool[1;;],
                Bool[1 1;],
                Bool[
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 1 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                ],
                Bool[
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                    0 0 1 1 0 0
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                ]
            )
        end
        @testset "1x{2,3} example" begin
            test_padding_for_unlike_images(
                Bool[1 1;],
                Bool[1 1 1;],
                Bool[
                    0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0
                    0 0 0 1 1 0 0 0
                    0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0
                ],
                Bool[
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 1 1 1 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                ]
            )
        end
        @testset "{2,3}x1 example" begin
            test_padding_for_unlike_images(
                Bool[
                    1;
                    1;;
                ],
                Bool[
                    1;
                    1;
                    1;;
                ],
                Bool[
                    0 0 0 0 0 0 0
                    0 0 0 0 0 0 0
                    0 0 0 0 0 0 0
                    0 0 0 1 0 0 0
                    0 0 0 1 0 0 0
                    0 0 0 0 0 0 0
                    0 0 0 0 0 0 0
                    0 0 0 0 0 0 0
                ],
                Bool[
                    0 0 0 0 0 0 0
                    0 0 0 0 0 0 0
                    0 0 0 0 0 0 0
                    0 0 0 1 0 0 0
                    0 0 0 1 0 0 0
                    0 0 0 1 0 0 0
                    0 0 0 0 0 0 0
                    0 0 0 0 0 0 0
                    0 0 0 0 0 0 0
                ]
            )
        end
    end


end