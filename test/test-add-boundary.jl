@testitem "add_boundary! basic functionality" begin
    using IceFloeTracker.Tracking: add_boundary!
    using IceFloeTracker.Segmentation: regionprops_table
    using DataFrames: DataFrame, nrow

    # Create a simple labeled image with two floes
    img = zeros(Int, 10, 10)
    img[2:4, 2:4] .= 1  # First floe (3x3)
    img[6:8, 6:8] .= 2  # Second floe (3x3)

    # Get regionprops with masks
    props = regionprops_table(img, properties=[:label, :area, :mask])

    # Add boundaries
    add_boundary!(props)

    # Test 1: Boundary column exists
    @test "boundary" ∈ names(props)

    # Test 2: Each boundary is a Matrix{Float64}
    for boundary in props.boundary
        @test boundary isa Matrix{Float64}
        @test size(boundary, 2) == 2  # [x y] format
        @test size(boundary, 1) >= 2   # At least 2 points
    end

    # Test 3: Number of boundaries matches number of floes
    @test nrow(props) == 2
end

@testitem "add_boundary! default reduc_factor" begin
    using IceFloeTracker.Tracking: add_boundary!
    using IceFloeTracker.Segmentation: regionprops_table
    using IceFloeTracker.Tracking: bwtraceboundary, resample_boundary

    # Create a labeled image
    img = zeros(Int, 12, 12)
    img[2:10, 2:10] .= 1  # 9x9 floe

    props = regionprops_table(img, properties=[:label, :mask])
    add_boundary!(props)  # Default reduc_factor=2

    # Compare with manual computation
    mask = props.mask[1]
    bd_traced = bwtraceboundary(mask)
    # Handle both Vector{CartesianIndex} and Vector{Vector{CartesianIndex}} returns
    bd_single = isa(bd_traced, Vector{Vector{CartesianIndex}}) ? bd_traced[1] : bd_traced
    bd_resampled = resample_boundary(bd_single, 2)  # reduc_factor=2

    # Boundaries should have similar number of points
    @test size(props.boundary[1], 1) == size(bd_resampled, 1)
end

@testitem "add_boundary! custom reduc_factor" begin
    using IceFloeTracker.Tracking: add_boundary!
    using IceFloeTracker.Segmentation: regionprops_table

    img = zeros(Int, 12, 12)
    img[2:10, 2:10] .= 1

    props1 = regionprops_table(img, properties=[:label, :mask])
    add_boundary!(props1; reduc_factor=1)  # No reduction

    props2 = regionprops_table(img, properties=[:label, :mask])
    add_boundary!(props2; reduc_factor=4)  # More reduction

    # With reduc_factor=1, should have more points than reduc_factor=4
    n_points_1 = size(props1.boundary[1], 1)
    n_points_4 = size(props2.boundary[1], 1)

    @test n_points_1 >= n_points_4
end

@testitem "add_boundary! boundary format" begin
    using IceFloeTracker.Tracking: add_boundary!
    using IceFloeTracker.Segmentation: regionprops_table

    # Create a simple square
    img = zeros(Int, 8, 8)
    img[2:6, 2:6] .= 1  # 5x5 square

    props = regionprops_table(img, properties=[:label, :mask])
    add_boundary!(props)

    boundary = props.boundary[1]

    # Test format
    @test eltype(boundary) == Float64
    @test size(boundary, 2) == 2

    # First and last points should be close (closed boundary)
    @test isapprox(boundary[1, :], boundary[end, :]; atol=1e-10)
end

@testitem "add_boundary! multiple floes" begin
    using IceFloeTracker.Tracking: add_boundary!
    using IceFloeTracker.Segmentation: regionprops_table, nrow

    # Create image with 3 floes of different sizes
    img = zeros(Int, 15, 15)
    img[2:4, 2:4] .= 1      # Small (3x3)
    img[6:10, 6:10] .= 2    # Medium (5x5)
    img[11:14, 11:14] .= 3  # Large (4x4)

    props = regionprops_table(img, properties=[:label, :area, :mask])
    add_boundary!(props)

    # All 3 floes should have boundaries
    @test nrow(props) == 3
    @test all(typeof.(props.boundary) .== Matrix{Float64})

    # Larger floes should generally have more boundary points
    n_points = size.(props.boundary, 1)
    @test n_points[2] >= n_points[1]  # Medium >= Small
    @test n_points[3] >= n_points[1]  # Large >= Small
end

@testitem "add_boundary! with add_floemasks workflow" begin
    using IceFloeTracker.Tracking: add_boundary!, add_floemasks!
    using IceFloeTracker.Segmentation: regionprops_table

    # Create labeled image
    img = zeros(Int, 10, 10)
    img[2:4, 2:4] .= 1
    img[6:8, 6:8] .= 2

    # Workflow: regionprops → add masks → add boundaries
    props = regionprops_table(img, properties=[:label, :area])
    add_floemasks!(props, img)
    add_boundary!(props)

    # Should have both masks and boundaries
    @test "mask" ∈ names(props)
    @test "boundary" ∈ names(props)
    @test nrow(props) == 2
    @test all(typeof.(props.boundary) .== Matrix{Float64})
end

@testitem "add_boundary! preserves other columns" begin
    using IceFloeTracker.Tracking: add_boundary!
    using IceFloeTracker.Segmentation: regionprops_table

    img = zeros(Int, 10, 10)
    img[2:6, 2:6] .= 1

    props = regionprops_table(img, properties=[:label, :area, :perimeter, :mask])
    original_cols = names(props)

    add_boundary!(props)

    # Original columns should still be there
    for col in original_cols
        @test col ∈ names(props)
    end

    # Plus the new boundary column
    @test "boundary" ∈ names(props)
end

@testitem "regionprops :mask deprecation warning" begin
    using IceFloeTracker.Segmentation: regionprops_table
    using Test

    img = zeros(Int, 10, 10)
    img[2:6, 2:6] .= 1

    # Verify that requesting :mask emits a deprecation warning
    @test_logs (:warn, r":mask.*deprecated") regionprops_table(img, properties=[:label, :mask])
end
