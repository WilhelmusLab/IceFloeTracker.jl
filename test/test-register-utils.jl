
@testitem "registration utilities" begin
    using IceFloeTracker: align_centroids, compute_centroid, shape_difference_rotation

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

        @testset "larger images" begin
            @test compute_centroid(
                Bool[
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
                ];
                rounded=false,
            ) == (9.0, 9.0)
            @test compute_centroid(
                Bool[
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
                ];
                rounded=false,
            ) == (8.964285714285714, 6.375)
        end
    end

    @testset "align centroids" begin
        @testset "simple cases" begin
            @test align_centroids(
                Bool[
                    1 0
                    0 0
                ],
                Bool[
                    0 0
                    0 1
                ],
            ) == (
                Bool[
                    0 0 0
                    0 1 0
                    0 0 0
                ],
                Bool[
                    0 0 0
                    0 1 0
                    0 0 0
                ],
            )
            @test align_centroids(
                Bool[1;;],
                Bool[
                    0 0
                    0 1
                ],
            ) == (
                Bool[
                    0 0
                    0 1
                ],
                Bool[
                    0 0
                    0 1
                ],
            )
            @test align_centroids(
                Bool[1 1;],
                Bool[
                    0 0
                    0 1
                ],
            ) == (
                Bool[
                    0 0 0
                    0 1 1
                ],
                Bool[
                    0 0 0
                    0 1 0
                ],
            )
            @test align_centroids(
                Bool[
                    0 0 0 0 0 1 1 1;
                ],
                Bool[
                    1 0 0;
                ],
            ) == (
                Bool[
                    0 0 0 0 0 1 1 1 0;
                ],
                Bool[
                    0 0 0 0 0 0 1 0 0;
                ],
            )
        end
        @testset "already centered image stays unchanged" begin
            function test_same_already_centered_image_is_unchanged(image)
                actual_result = align_centroids(image, image)
                expected_result = (image, image)
                @test actual_result == expected_result
            end

            test_same_already_centered_image_is_unchanged(Bool[1;;])
            test_same_already_centered_image_is_unchanged(Bool[1 1;])
            test_same_already_centered_image_is_unchanged(Bool[0 1 0;])
            test_same_already_centered_image_is_unchanged(Bool[1 1 1 1;])
            test_same_already_centered_image_is_unchanged(Bool[
                1; 1;;
            ])
            test_same_already_centered_image_is_unchanged(Bool[
                0; 1; 0;;
            ])
            test_same_already_centered_image_is_unchanged(Bool[0 1 1 1 0;])
            test_same_already_centered_image_is_unchanged(Bool[0 0 1 1 1 0 0;])
            test_same_already_centered_image_is_unchanged(Bool[
                0 0 1 1 1 0 0
                0 0 1 1 1 0 0
            ])
            test_same_already_centered_image_is_unchanged(
                Bool[
                    0 0 1 1 1 0 0
                    0 0 1 1 1 0 0
                    0 0 1 1 1 0 0
                ]
            )
            test_same_already_centered_image_is_unchanged(
                Bool[
                    0 0 0 0 0 0 0
                    0 0 1 1 1 0 0
                    0 0 1 1 1 0 0
                    0 0 1 1 1 0 0
                    0 0 0 0 0 0 0
                ]
            )
            test_same_already_centered_image_is_unchanged(
                Bool[
                    0 0 0 0 0 0 0
                    0 0 0 0 0 0 0
                    0 0 1 1 1 0 0
                    0 0 1 1 1 0 0
                    0 0 1 1 1 0 0
                    0 0 0 0 0 0 0
                    0 0 0 0 0 0 0
                ],
            )
            test_same_already_centered_image_is_unchanged(
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
                ],
            )
        end
        @testset "offset image" begin
            @test align_centroids(
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
                ],
                Bool[
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 1 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                ],
            ) == (
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
                ],
                Bool[
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 1 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                ],
            )

            @test align_centroids(
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
                ],
                Bool[
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 1 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 1 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                ],
            ) == (
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
                ],
                Bool[
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 1 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 1 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0 0
                ],
            )
            @test align_centroids(
                Bool[
                    0 0 0 0 0
                    0 1 1 1 0
                    0 1 1 1 0
                    0 1 1 1 0
                    0 0 0 0 0
                ],
                Bool[
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 1 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                ],
            ) == (
                Bool[
                    0 0 0 0 0
                    0 1 1 1 0
                    0 1 1 1 0
                    0 1 1 1 0
                    0 0 0 0 0
                ],
                Bool[
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 1 0 0
                    0 0 0 0 0
                    0 0 0 0 0
                ],
            )
            @test align_centroids(
                Bool[
                    0 0 0 0 0
                    0 1 1 1 0
                    0 1 1 0 0
                    0 1 0 0 0
                    0 0 0 0 0
                ],
                Bool[
                    0 0 0 0 0
                    0 0 0 0 0
                    0 0 1 1 0
                    0 0 1 1 0
                    0 0 0 0 0
                ],
            ) == (
                Bool[
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                    0 0 1 1 1 0
                    0 0 1 1 0 0
                    0 0 1 0 0 0
                    0 0 0 0 0 0
                ],
                Bool[
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                    0 0 1 1 0 0
                    0 0 1 1 0 0
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                ],
            )
        end
        @testset "larger realistic image" begin
            im1 = Bool[
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 1 1 1 0 0 0
                0 0 0 0 0 0 0 0 1 1 1 1 0 0 0
                0 0 0 0 0 0 1 1 1 1 1 1 1 0 0
                0 0 0 0 1 1 1 1 1 1 1 1 1 0 0
                0 0 0 1 1 1 1 1 1 1 1 1 0 0 0
                0 0 1 1 1 1 1 1 1 1 1 0 0 0 0
                0 0 1 1 1 1 1 1 1 0 0 0 0 0 0
                0 0 0 1 1 1 1 0 0 0 0 0 0 0 0
                0 0 0 1 1 1 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            ]
            im2 = Bool[
                0 0 0 0 0 0 0 0 0 0 0 0 0
                0 1 0 0 1 0 0 1 1 0 1 1 0
                0 0 1 1 1 1 1 1 1 1 1 1 0
                0 1 1 1 1 1 1 1 1 1 1 1 0
                0 1 1 1 1 1 1 1 1 1 1 1 0
                0 0 1 1 1 1 1 1 1 1 1 1 0
                0 1 0 0 1 0 0 1 1 0 1 1 0
                0 0 0 0 0 0 0 0 0 0 0 0 0
            ]
            result = align_centroids(im1, im2)
            @test size(result[1]) == size(result[2])
        end
    end

    @testset "shape difference" begin
        no_rotation = [0.0]
        @test shape_difference_rotation(Bool[1;;], Bool[1;;], no_rotation) ==
            [(; angle=0.0, shape_difference=0)]
        @test shape_difference_rotation(
            Bool[
                1 1;
            ],
            Bool[1;;],
            no_rotation,
        ) == [(; angle=0.0, shape_difference=1)]
        @test shape_difference_rotation(
            Bool[
                1 1
                1 0
            ],
            Bool[1;;],
            no_rotation,
        ) == [(; angle=0.0, shape_difference=2)]

        @test shape_difference_rotation(
            Bool[
                1 1
                1 0
            ],
            Bool[
                1 1
                1 0
            ],
            no_rotation,
        ) == [(; angle=0.0, shape_difference=0)]
        @test shape_difference_rotation(
            Bool[
                1 1
                1 0
            ],
            Bool[
                1 1
                1 1
            ],
            no_rotation,
        ) == [(; angle=0.0, shape_difference=1)]
        @test shape_difference_rotation(
            Bool[
                0 1
                1 1
            ],
            Bool[
                1 1
                1 0
            ],
            [0.0, π],
        ) == [(; angle=0.0, shape_difference=2), (; angle=Float64(π), shape_difference=0)]
        @test shape_difference_rotation(
            Bool[
                1 1 1
                1 1 0
                1 0 0
            ],
            Bool[
                1 1
                1 1
            ],
            no_rotation,
        ) == [(; angle=0.0, shape_difference=2)]
        @test shape_difference_rotation(
            Bool[
                1 1 1
                1 1 0
                1 0 0
            ],
            Bool[1;;],
            no_rotation,
        ) == [(; angle=0.0, shape_difference=5)]
        @test shape_difference_rotation(
            Bool[
                0 0 0 0 0
                0 1 1 1 0
                0 1 1 0 0
                0 1 0 0 0
                0 0 0 0 0
            ],
            Bool[
                0 0 0 0 0
                0 0 0 0 0
                0 0 1 1 0
                0 0 1 1 0
                0 0 0 0 0
            ],
            no_rotation,
        ) == [(; angle=0.0, shape_difference=2)]
        @test shape_difference_rotation(
            Bool[
                1 1 1
                1 0 1
                1 1 1
            ],
            Bool[
                0 0 0
                0 1 0
                0 0 0
            ],
            no_rotation,
        ) == [(; angle=0.0, shape_difference=9)]
        @test shape_difference_rotation(
            Bool[
                0 0 0 0 0
                0 1 1 1 0
                0 1 0 1 0
                0 1 1 1 0
                0 0 0 0 0
            ],
            Bool[
                0 0 0 0 0
                0 0 0 0 0
                0 0 1 0 0
                0 0 0 0 0
                0 0 0 0 0
            ],
            no_rotation,
        ) == [(; angle=0.0, shape_difference=9)]

        @test shape_difference_rotation(
            Bool[
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0
                0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0
                0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0
                0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0
                0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0
                0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0
                0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0
                0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0
                0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            ],
            Bool[
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0
                0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0
                0 0 0 1 1 1 1 1 1 1 0 0 0 0 0 0
                0 0 0 0 1 1 1 1 1 1 1 0 0 0 0 0
                0 0 0 0 0 1 1 1 1 1 1 1 0 0 0 0
                0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 0
                0 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0
                0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0
                0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
            ],
            [0.0, π / 3.0, π / 2.0, 3π / 2.0, π],
        ) == [
            (; angle=0.0, shape_difference=58.0),
            (; angle=Float64(π / 3.0), shape_difference=18.0),
            (; angle=Float64(π / 2.0), shape_difference=14.0),
            (; angle=Float64(3π / 2.0), shape_difference=14.0),
            (; angle=Float64(π), shape_difference=58.0),
        ]
    end
end

@testsnippet CoarseToFineSetup begin
    using IceFloeTracker.Tracking:
        register,
        register_coarse_to_fine,
        register_coarse_angles_rad,
        register_default_angles_rad,
        imrotate_bin_clockwise_radians

    # an ellipse with an offset notch, padded so that no rotation used here clips it;
    # 1° rotations of this floe are distinguishable after imrotate_bin
    floe = falses(101, 101)
    for i in 1:101, j in 1:101
        x, y = j - 51, i - 51
        floe[i, j] = (x / 38)^2 + (y / 20)^2 <= 1 && !(x > 15 && y < -2)
    end
    rotated(deg) = imrotate_bin_clockwise_radians(floe, deg2rad(deg))
end

@testitem "register_coarse_angles_rad spans ±30° in 5° steps" setup = [CoarseToFineSetup] begin
    @test length(register_coarse_angles_rad) == 13
    @test minimum(register_coarse_angles_rad) ≈ -π / 6
    @test maximum(register_coarse_angles_rad) ≈ π / 6
    @test all(diff(sort(register_coarse_angles_rad)) .≈ π / 36)
    # tie order matches register_default_angles_rad: |angle| ascending, positive first
    @test issorted(register_coarse_angles_rad; by=θ -> (abs(θ), -θ))
end

@testitem "register_coarse_to_fine recovers rotations to 1°" setup = [CoarseToFineSetup] begin
    for deg in (-13, -7, -3, 0, 2, 4, 8, 13, 22)
        @test register_coarse_to_fine(floe, rotated(deg)) ≈ deg2rad(deg) atol = 1e-9
    end
    # the 5° full grid cannot resolve the same rotation
    @test !isapprox(register(floe, rotated(13)), deg2rad(13); atol=deg2rad(1))
end

@testitem "register_coarse_to_fine rotates fewer times than the full grid" setup = [
    CoarseToFineSetup
] begin
    n = Ref(0)
    counting_rotate(im, θ) = (n[] += 1; imrotate_bin_clockwise_radians(im, θ))
    register_coarse_to_fine(floe, rotated(8); imrotate_function=counting_rotate)
    @test n[] <= 22
    @test n[] < length(register_default_angles_rad)
end

@testitem "register_coarse_to_fine breaks ties toward 0" setup = [CoarseToFineSetup] begin
    # every rotation of a centered disc scores the same
    disc = [hypot(i - 26, j - 26) <= 20 for i in 1:51, j in 1:51]
    @test register_coarse_to_fine(disc, disc) == 0
end

@testitem "register_coarse_to_fine warns when the coarse optimum is on the window edge" setup = [
    CoarseToFineSetup
] begin
    @test_logs (:warn, r"edge") register_coarse_to_fine(floe, rotated(60))
    @test_logs register_coarse_to_fine(floe, rotated(13))
end
