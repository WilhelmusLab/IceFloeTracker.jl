using IceFloeTracker: pad_images, compute_centroid

@testset "registration utilities" begin
    @testset "centroid" begin
        @testset "zeroes" begin
            @test isequal(compute_centroid(Bool[0;;]; rounded=false), (NaN, NaN))
            @test_throws InexactError compute_centroid(Bool[0;;]; rounded=true)
        end
        @testset "one row" begin
            @test compute_centroid(Bool[1;;]; rounded=false) == (1.0, 1.0)
            @test compute_centroid(Bool[1;;]; rounded=true) == (1, 1)
            @test compute_centroid(Bool[1 1;]; rounded=false) == (1.0, 1.5)
            @test compute_centroid(Bool[1 1;]; rounded=true) == (1, 2)
            @test compute_centroid(Bool[1 1 1;]; rounded=false) == (1.0, 2.0)
            @test compute_centroid(Bool[1 1 1;]; rounded=true) == (1, 2)
            @test compute_centroid(Bool[1 1 1 1;]; rounded=false) == (1.0, 2.5)
            @test compute_centroid(Bool[1 1 1 1;]; rounded=true) == (1, 2)  # this test seems a bit flaky
            @test compute_centroid(Bool[1 1 1 1 1;]; rounded=false) == (1.0, 3.0)
            @test compute_centroid(Bool[1 1 1 1 1;]; rounded=true) == (1, 3)
        end

        @testset "one column" begin
            @test compute_centroid(Bool[1; 1;;]; rounded=false) == (1.5, 1.0)
            @test compute_centroid(Bool[1; 1;;]; rounded=true) == (2, 1)
            @test compute_centroid(Bool[1; 1; 1;;]; rounded=false) == (2.0, 1.0)
            @test compute_centroid(Bool[1; 1; 1;;]; rounded=true) == (2, 1)
            @test compute_centroid(Bool[1; 1; 1; 1;;]; rounded=false) == (2.5, 1.0)
            @test compute_centroid(Bool[1; 1; 1; 1;;]; rounded=true) == (2, 1)
            @test compute_centroid(Bool[1; 1; 1; 1; 1;;]; rounded=false) == (3.0, 1.0)
            @test compute_centroid(Bool[1; 1; 1; 1; 1;;]; rounded=true) == (3, 1)
        end

        @testset "simple matrix" begin
            @test compute_centroid(Bool[
                0 0 0
                0 1 0
                0 0 0
            ]) == (2, 2)
            @test compute_centroid(Bool[
                1 0 1
                0 1 0
                1 0 1
            ]) == (2, 2)
            @test compute_centroid(Bool[
                1 1 1
                1 1 1
                1 1 1
            ]) == (2, 2)
            @test compute_centroid(Bool[
                1 0 0
                0 1 0
                0 0 1
            ]) == (2, 2)
            @test compute_centroid(Bool[
                0 0 1
                0 1 0
                1 0 0
            ]) == (2, 2)
            @test compute_centroid(Bool[
                1 0 1
                0 0 0
                1 0 1
            ]) == (2, 2)
            @test compute_centroid(Bool[
                1 1 1
                1 0 1
                1 1 1
            ]) == (2, 2)
            @test compute_centroid(Bool[
                0 1 0
                1 0 1
                0 1 0
            ]) == (2, 2)
            @test compute_centroid(Bool[
                1 0 0
                0 0 0
                0 0 1
            ]) == (2, 2)
        end

        @testset "scaled matrices" begin
            function test_measured_centroid_is_independent_of_scaling(array, scaling_factors; rounded)
                for f in scaling_factors
                    @test compute_centroid(array; rounded) == compute_centroid(array * f; rounded)
                end
            end
            test_measured_centroid_is_independent_of_scaling(
                abs.(randn(3, 3)),
                abs.(randn(100));
                rounded=true
            )
            test_measured_centroid_is_independent_of_scaling(
                abs.(randn(5, 5)),
                abs.(randn(100));
                rounded=false
            )
        end

        @testset "larger images" begin
            @test compute_centroid(Bool[
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0
                    0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0
                    0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0
                    0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0
                    0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0
                    0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0
                    0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0
                    0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0
                    0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0
                    0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0
                    0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                ]; rounded=false) == (9.0, 9.0)
            @test compute_centroid(Bool[
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0
                    0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0
                    0 0 0 0 1 0 0 0 1 1 1 1 0 0 0 0
                    0 0 0 1 1 1 0 1 1 1 1 0 0 0 0 0
                    0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0
                    1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
                    0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0
                    0 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0
                    0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0
                    0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0
                    0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                ]; rounded=false) == (8.964285714285714, 6.375)

        end
    end

    @testset "padding" begin
        @testset "two identical inputs are padded the same" begin
            function test_padding_for_identical_inputs(image, expected_result)
                result = pad_images(image, image)
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
                @test result1 == expected_result1
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
end