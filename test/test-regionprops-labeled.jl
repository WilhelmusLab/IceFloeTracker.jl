@testset "regionprops-labeled.jl" begin
    println("------------------------------------------------------------")
    println("-------------- regionprops (labeled) Tests ----------------")

    img1 = Int[
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   # row  1
        0 1 1 0 2 2 2 0 3 3 3 3 0 0 0 0   #      2
        0 1 1 0 2 2 2 0 3 3 3 3 0 0 0 0   #      3
        0 0 0 0 0 2 2 0 3 3 3 3 0 0 0 0   #      4
        0 0 0 0 0 0 0 0 3 3 3 0 0 0 0 0   #      5
        0 4 4 4 4 5 0 0 0 0 0 0 0 0 0 0   #      6
        0 4 5 5 5 5 0 0 0 0 0 0 0 0 0 0   #      7
        0 4 5 5 0 0 0 0 0 0 0 0 0 0 0 0   #      8
        0 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0   #      9
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   #     10
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   #     11
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   #     12
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   #     13
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   #     14
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   #     15
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   #     16
        #= col
        1 2 3 4 5 6 7 8 9 a b c d e f g
        =#

    ]

    @info "testing props with label"
    props_with_label = DataFrame([
        1 1 1 4 4
        2 1 4 5 8
        3 1 8 6 13
        4 5 1 10 7
        5 5 1 10 7
    ], [:label, :min_row, :min_col, :max_row, :max_col])


    @test cropfloe(img1, props_with_label, 1) == [
        0 0 0 0
        0 1 1 0
        0 1 1 0
        0 0 0 0
    ]

    @test cropfloe(img1, props_with_label, 2) == [
        0 0 0 0 0
        0 1 1 1 0
        0 1 1 1 0
        0 0 1 1 0
        0 0 0 0 0
    ]

    @test cropfloe(img1, props_with_label, 3) == [
        0 0 0 0 0 0
        0 1 1 1 1 0
        0 1 1 1 1 0
        0 1 1 1 1 0
        0 1 1 1 0 0
        0 0 0 0 0 0
    ]

    @test cropfloe(img1, props_with_label, 4) == [
        0 0 0 0 0 0 0
        0 1 1 1 1 0 0
        0 1 0 0 0 0 0
        0 1 0 0 0 0 0
        0 1 1 1 1 1 0
        0 0 0 0 0 0 0
    ]

    @test cropfloe(img1, props_with_label, 5) == [
        0 0 0 0 0 0 0
        0 0 0 0 0 1 0
        0 0 1 1 1 1 0
        0 0 1 1 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
    ]

    @info "testing props without label"
    props_without_label = DataFrame([
        1 1 4 4
        1 4 5 8
        1 8 6 13
        5 1 10 7
        5 1 10 7
    ], [:min_row, :min_col, :max_row, :max_col])

    @test cropfloe(img1, props_without_label, 1) == [
        0 0 0 0
        0 1 1 0
        0 1 1 0
        0 0 0 0
    ]

    @test cropfloe(img1, props_without_label, 2) == [
        0 0 0 0 0
        0 1 1 1 0
        0 1 1 1 0
        0 0 1 1 0
        0 0 0 0 0
    ]

    @test cropfloe(img1, props_without_label, 3) == [
        0 0 0 0 0 0
        0 1 1 1 1 0
        0 1 1 1 1 0
        0 1 1 1 1 0
        0 1 1 1 0 0
        0 0 0 0 0 0
    ]

    # This can't distinguish between the two interlocking regions, so it merges them
    @test cropfloe(img1, props_without_label, 4) == [
        0 0 0 0 0 0 0
        0 1 1 1 1 1 0
        0 1 1 1 1 1 0
        0 1 1 1 0 0 0
        0 1 1 1 1 1 0
        0 0 0 0 0 0 0
    ]

    @test cropfloe(img1, props_without_label, 4) == cropfloe(img1, props_without_label, 5)



    @info "testing props with only label"
    props_only_label = DataFrame(label=1:5)

    @test cropfloe(img1, props_only_label, 1) == [
        0 0 0 0
        0 1 1 0
        0 1 1 0
        0 0 0 0
    ]

    @test cropfloe(img1, props_only_label, 2) == [
        0 0 0 0 0
        0 1 1 1 0
        0 1 1 1 0
        0 0 1 1 0
        0 0 0 0 0
    ]

    @test cropfloe(img1, props_only_label, 3) == [
        0 0 0 0 0 0
        0 1 1 1 1 0
        0 1 1 1 1 0
        0 1 1 1 1 0
        0 1 1 1 0 0
        0 0 0 0 0 0
    ]

    @test cropfloe(img1, props_only_label, 4) == [
        0 0 0 0 0 0 0
        0 1 1 1 1 0 0
        0 1 0 0 0 0 0
        0 1 0 0 0 0 0
        0 1 1 1 1 1 0
        0 0 0 0 0 0 0
    ]

    @test cropfloe(img1, props_only_label, 5) == [
        0 0 0 0 0 0
        0 0 0 0 1 0
        0 1 1 1 1 0
        0 1 1 0 0 0
        0 0 0 0 0 0
    ]


    @info "testing values with label"
    # cropfloe
    @test cropfloe(img1, 1, 1, 4, 4, 1) == [
        0 0 0 0
        0 1 1 0
        0 1 1 0
        0 0 0 0
    ]

    @test cropfloe(img1, 1, 4, 5, 8, 2) == [
        0 0 0 0 0
        0 1 1 1 0
        0 1 1 1 0
        0 0 1 1 0
        0 0 0 0 0
    ]

    @test cropfloe(img1, 1, 8, 6, 13, 3) == [
        0 0 0 0 0 0
        0 1 1 1 1 0
        0 1 1 1 1 0
        0 1 1 1 1 0
        0 1 1 1 0 0
        0 0 0 0 0 0
    ]

    @test cropfloe(img1, 5, 1, 10, 7, 4) == [
        0 0 0 0 0 0 0
        0 1 1 1 1 0 0
        0 1 0 0 0 0 0
        0 1 0 0 0 0 0
        0 1 1 1 1 1 0
        0 0 0 0 0 0 0
    ]

    @test cropfloe(img1, 5, 1, 10, 7, 5) == [
        0 0 0 0 0 0 0
        0 0 0 0 0 1 0
        0 0 1 1 1 1 0
        0 0 1 1 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
    ]

    @test cropfloe(img1, 1) == [
        0 0 0 0
        0 1 1 0
        0 1 1 0
        0 0 0 0
    ]

    @test cropfloe(img1, 2) == [
        0 0 0 0 0
        0 1 1 1 0
        0 1 1 1 0
        0 0 1 1 0
        0 0 0 0 0
    ]

    @test cropfloe(img1, 3) == [
        0 0 0 0 0 0
        0 1 1 1 1 0
        0 1 1 1 1 0
        0 1 1 1 1 0
        0 1 1 1 0 0
        0 0 0 0 0 0
    ]

    @test cropfloe(img1, 4) == [
        0 0 0 0 0 0 0
        0 1 1 1 1 0 0
        0 1 0 0 0 0 0
        0 1 0 0 0 0 0
        0 1 1 1 1 1 0
        0 0 0 0 0 0 0
    ]

    @test cropfloe(img1, 5) == [
        0 0 0 0 0 0
        0 0 0 0 1 0
        0 1 1 1 1 0
        0 1 1 0 0 0
        0 0 0 0 0 0
    ]

end
