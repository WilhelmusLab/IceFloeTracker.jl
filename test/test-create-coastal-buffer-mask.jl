
@testitem "Create Coastal Buffer Mask" begin
    using Images: centered

    function coastal_buffer_creation(landmask, struct_elem, expected_mask; kwargs...)
        actual_mask = create_coastal_buffer_mask(landmask, struct_elem; kwargs...)
        if actual_mask != expected_mask
            @info "expected_mask:"
            display(expected_mask)
            @info "actual_mask:"
            display(actual_mask)
            @info "difference:"
            display(expected_mask .!= actual_mask)
        end
        return actual_mask == expected_mask
    end

    @test coastal_buffer_creation(
        BitMatrix([
            0 0 0 0 0
            0 1 1 1 0
            0 1 0 1 0
            0 1 1 1 0
            0 0 0 0 0
        ]),
        trues(3, 3),
        BitMatrix([
            1 1 1 1 1
            1 1 1 1 1
            1 1 1 1 1
            1 1 1 1 1
            1 1 1 1 1
        ]);
        fill_holes_min_pixels=0,
        fill_holes_max_pixels=0,
    )

    @test coastal_buffer_creation(
        BitMatrix([
            0 0 1 0 0
            0 0 1 0 0
            0 0 1 0 0
            0 0 1 0 0
            0 0 1 0 0
        ]),
        trues(3, 3),
        BitMatrix([
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 1 0
        ]);
        fill_holes_min_pixels=0,
        fill_holes_max_pixels=0,
    )

    @test coastal_buffer_creation(
        BitMatrix([
            0 0 0 0 0 0 0
            0 1 1 1 0 0 0
            0 1 0 1 1 0 0
            0 1 1 1 0 0 0
            0 0 0 0 0 0 0
        ]),
        trues(3, 3),
        BitMatrix([
            1 1 1 1 1 0 0
            1 1 1 1 1 1 0
            1 1 1 1 1 1 0
            1 1 1 1 1 1 0
            1 1 1 1 1 0 0
        ]);
        fill_holes_min_pixels=0,
        fill_holes_max_pixels=0,
    )
    @test coastal_buffer_creation(
        BitMatrix(
            [
                0 0 0 0 0 0 0 0
                0 0 1 0 1 0 0 0
                0 0 1 1 0 1 0 0
                0 0 1 0 0 0 1 0
                0 0 1 0 0 0 1 0
                0 0 1 0 1 0 0 0
                0 0 0 0 0 0 0 0
            ],
        ),
        trues(3, 3),
        [
            0 1 1 1 1 1 0 0
            0 1 1 1 1 1 1 0
            0 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1
            0 1 1 1 1 1 1 1
            0 1 1 1 1 1 0 0
        ];
        fill_holes_min_pixels=0,
        fill_holes_max_pixels=0,
    )

    @test coastal_buffer_creation(
        BitMatrix(
            [
                0 0 0 0 0 0 0 0
                0 0 1 0 1 0 0 0
                0 0 1 1 0 1 0 0
                0 0 1 0 0 0 1 0
                0 0 1 0 0 0 1 0
                0 0 1 0 1 0 0 0
                0 0 0 0 0 0 0 0
            ],
        ),
        BitMatrix([
            0 1 0
            0 1 0
            0 1 0
        ]),
        BitMatrix(
            [
                0 0 1 1 1 0 0 0
                0 0 1 1 1 1 0 0
                0 0 1 1 1 1 1 0
                0 0 1 1 1 1 1 0
                0 0 1 0 1 0 1 0
                0 0 1 0 1 0 1 0
                0 0 1 0 1 0 0 0
            ],
        );
        fill_holes_min_pixels=0,
        fill_holes_max_pixels=1,
    )

    @test coastal_buffer_creation(
        BitMatrix([
            1 1 1 1 1
            1 0 0 0 1
            1 0 0 0 1
            1 0 0 0 1
            1 1 1 1 1
        ]),
        trues(1, 1),
        BitMatrix([
            1 1 1 1 1
            1 1 1 1 1
            1 1 1 1 1
            1 1 1 1 1
            1 1 1 1 1
        ]);
        fill_holes_min_pixels=0,
        fill_holes_max_pixels=25,
    )

    @test coastal_buffer_creation(
        BitMatrix([
            1 1 1 1 1
            1 0 0 0 1
            1 0 0 0 1
            1 0 0 0 1
            1 1 1 1 1
        ]),
        trues(1, 1),
        BitMatrix([
            1 1 1 1 1
            1 0 0 0 1
            1 0 0 0 1
            1 0 0 0 1
            1 1 1 1 1
        ]);
        fill_holes_min_pixels=0,
        fill_holes_max_pixels=8,
    )
end