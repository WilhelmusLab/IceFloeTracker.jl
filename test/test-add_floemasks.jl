@testitem "add_floemasks!" begin
    import Images: SegmentedImage, labels_map
    @testset "SegmentedImage" begin
        labels = [
            0 0 0 0 0 0 0
            0 6 6 0 0 2 0
            0 6 6 0 0 2 2
            0 0 0 0 0 0 0
            0 3 3 0 4 4 4
            0 0 3 0 4 4 4
        ]
        seg_img = SegmentedImage(labels, labels)
        props = regionprops_table(seg_img; properties=[:label, :mask])
        # add_floemasks!(props, seg_img)

        @test props[1, :].label == 2
        @test props[1, :].mask == BitMatrix([
            1 0
            1 1
        ])

        @test props[2, :].label == 3
        @test props[2, :].mask == BitMatrix([
            1 1
            0 1
        ])

        @test props[3, :].label == 4
        @test props[3, :].mask == BitMatrix([
            1 1 1
            1 1 1
        ])

        @test props[4, :].label == 6
        @test props[4, :].mask == BitMatrix([
            1 1
            1 1
        ])
    end

    @testset "Vector{SegmentedImage}" begin
        labels::Vector{Matrix{Int}} = [
            [
                0 0 0 0 0 0 0
                0 6 6 0 0 2 0
                0 6 6 0 0 2 2
                0 0 0 0 0 0 0
                0 3 3 0 4 4 4
                0 0 3 0 4 4 4
            ],
            [
                0 0 0 0 0 0 0
                0 5 5 0 3 0 0
                0 5 5 0 3 3 0
                0 0 0 0 0 0 0
                0 2 2 4 4 4 0
                0 0 2 4 4 4 0
            ],
        ]
        seg_imgs = SegmentedImage.(labels, labels)
        props = regionprops_table.(seg_imgs; properties=[:label, :mask])
        # add_floemasks!.(props, seg_imgs)


        @test props[1][1, :].label == 2
        @test props[1][1, :].mask == BitMatrix([
            1 0
            1 1
        ])

        @test props[1][2, :].label == 3
        @test props[1][2, :].mask == BitMatrix([
            1 1
            0 1
        ])

        @test props[1][3, :].label == 4
        @test props[1][3, :].mask == BitMatrix([
            1 1 1
            1 1 1
        ])

        @test props[1][4, :].label == 6
        @test props[1][4, :].mask == BitMatrix([
            1 1
            1 1
        ])

        @test props[2][1, :].label == 2
        @test props[2][1, :].mask == BitMatrix([
            1 1
            0 1
        ])

        @test props[2][2, :].label == 3
        @test props[2][2, :].mask == BitMatrix([
            1 0
            1 1
        ])

        @test props[2][3, :].label == 4
        @test props[2][3, :].mask == BitMatrix([
            1 1 1
            1 1 1
        ])

        @test props[2][4, :].label == 5
        @test props[2][4, :].mask == BitMatrix([
            1 1
            1 1
        ])

    end
end