@testsnippet RegionpropsFixture begin
    using Random
    import Images: label_components

    Random.seed!(123)
    bw_img = Bool.(rand([0, 1], 5, 10))
    bw_img[end, 7] = 1
    label_img = label_components(bw_img, trues(3, 3))
end

@testitem "regionprops" setup = [RegionpropsFixture] begin
    import DataFrames: DataFrame, nrow
    import IceFloeTracker: PixelConvexArea, PolygonConvexArea

    properties = [
        "label",
        "centroid",
        "area",
        "major_axis_length",
        "minor_axis_length",
        "bbox",
        "perimeter",
        "orientation",
        "circularity",
        "solidity",
    ]
    extra_props = nothing

    table = regionprops_table(label_img, bw_img; properties=properties, minimum_area=1)
    total_labels = maximum(label_img)

    # Tests for regionprops_table
    @test typeof(table) <: DataFrame # check correct data type
    @test nrow(table) == total_labels - 1 # One point has area 1, we exclude that point intentionally
    @test all(
        [
            "area",
            "min_row",
            "min_col",
            "max_row",
            "max_col",
            "row_centroid",
            "col_centroid",
            "label",
            "major_axis_length",
            "minor_axis_length",
            "orientation",
            "perimeter",
            "circularity",
            "solidity",
        ] ⊆ names(table),
    )
    # check all requested properties are present

    # Tests for regionprops
    regions = regionprops(label_img, bw_img)
    @test table.area[2] == regions[:area][2]

    # Check floe masks generation and correct cropping
    add_floemasks!(table, label_img)
    @test all(
        [
            length(unique(label_components(table.mask[i], trues(3, 3)))) for
            i in 1:nrow(table)
        ] .== [2, 2],
    )

    # Test that algorithm options at least run
    regionprops(label_img; properties=["perimeter"], perimeter_algorithm=BenkridCrookes())
    regionprops(
        label_img; properties=["convex_area"], convex_area_algorithm=PolygonConvexArea()
    )
    regionprops(
        label_img; properties=["convex_area"], convex_area_algorithm=PixelConvexArea()
    )

    # Test that we can get complex measures without listing everything
    regionprops(
        label_img; properties=["area", "perimeter", "convex_area", "circularity", "solidity"], convex_area_algorithm=PixelConvexArea()
    )

    # Regression test: cross-shaped region with 4 pixels should not error
    # (see https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/919)
    cross_img = [
        0 0 0 0
        0 1 0 0
        0 1 1 0
        0 1 0 0
        0 0 0 0
    ]
    @test_nowarn regionprops_table(
        cross_img; properties=["convex_area"], convex_area_algorithm=PixelConvexArea()
    )
    @test_nowarn regionprops_table(
        cross_img; properties=["convex_area"], convex_area_algorithm=PolygonConvexArea()
    )
    cross_result = regionprops_table(
        cross_img; properties=["convex_area"], convex_area_algorithm=PixelConvexArea()
    )
    @test nrow(cross_result) == 1

    @test isnan(cross_result.convex_area[1])
end

@testitem "regionprops: circularity and solidity work without their respective base properties" setup = [
    RegionpropsFixture,
] begin
    area = regionprops(label_img; properties=[:area])[:area]
    perimeter = regionprops(label_img; properties=[:perimeter])[:perimeter]
    convex_area = regionprops(label_img; properties=[:convex_area])[:convex_area]

    circularity_only = regionprops(label_img; properties=[:circularity])
    @test collect(keys(circularity_only)) == [:circularity] # :perimeter isn't leaked into the output
    @test circularity_only[:circularity] == 4 .* π .* area ./ perimeter .^ 2

    solidity_only = regionprops(label_img; properties=[:solidity])
    @test collect(keys(solidity_only)) == [:solidity] # :convex_area isn't leaked into the output
    @test solidity_only[:solidity] == area ./ convex_area
end

@testitem "regionprops table output should include all implied columns even if there are no rows" begin
    using DataFrames: nrow

    label_img = zeros(Int, 5, 5)
    @testset "defaults" begin
        table = regionprops_table(label_img)
        @test "min_row" ∈ names(table)
        @test "max_row" ∈ names(table)
        @test "min_col" ∈ names(table)
        @test "max_col" ∈ names(table)
        @test "row_centroid" ∈ names(table)
        @test "col_centroid" ∈ names(table)
        @test "label" ∈ names(table)
        @test "area" ∈ names(table)
        @test "major_axis_length" ∈ names(table)
        @test "minor_axis_length" ∈ names(table)
        @test "convex_area" ∈ names(table)
        @test "perimeter" ∈ names(table)
        @test "orientation" ∈ names(table)
    end

    @testset "bbox" begin
        table = regionprops_table(label_img; properties=["bbox"])
        @test "min_row" ∈ names(table)
        @test "max_row" ∈ names(table)
        @test "min_col" ∈ names(table)
        @test "max_col" ∈ names(table)
    end

    @testset "centroid" begin
        table = regionprops_table(label_img; properties=["centroid"])
        @test "row_centroid" ∈ names(table)
        @test "col_centroid" ∈ names(table)
    end
end

@testitem "addlatlon! adds latitude and longitude columns to regionprops table" begin
    using IceFloeTracker
    import DataFrames: DataFrame
    using FileIO

    refimage = "test_inputs/latlon/latlon_test_image-2020-06-21T00_00_00Z.tif"
    data = load(refimage)
    nrows, ncols = size(data)

    propdf = DataFrame(;
        row_centroid=[1.0, 1.0, Float64(nrows), Float64(nrows)],
        col_centroid=[1.0, Float64(ncols), 1.0, Float64(ncols)],
        area=ones(4),
        convex_area=ones(4),
        minor_axis_length=ones(4),
        major_axis_length=ones(4),
        perimeter=ones(4),
    )

    addlatlon!(propdf, refimage)

    @test all(["latitude", "longitude"] ⊆ names(propdf)) # check new columns added
    @info propdf
    top_left, top_right, bottom_left, bottom_right = 1, 2, 3, 4

    # For this georeferenced image, northern pixels (top rows) should have higher latitude.
    @test propdf.latitude[top_left] > propdf.latitude[bottom_left]
    @test propdf.latitude[top_right] > propdf.latitude[bottom_right]

    # Longitudes should increase from left to right at both top and bottom rows.
    @test propdf.longitude[top_right] > propdf.longitude[top_left]
    @test propdf.longitude[bottom_right] > propdf.longitude[bottom_left]
end

@testitem "converttounits maps corner centroids to sensible lat/lon ordering" begin
    using IceFloeTracker
    import DataFrames: DataFrame

    refimage = "test_inputs/latlon/latlon_test_image-2020-06-21T00_00_00Z.tif"
    latlondata = latlon(refimage)
    latitude = latlondata[:latitude]
    nrows, ncols = size(latitude)

    propdf = DataFrame(;
        row_centroid=[1.0, 1.0, Float64(nrows), Float64(nrows)],
        col_centroid=[1.0, Float64(ncols), 1.0, Float64(ncols)],
        area=ones(4),
        convex_area=ones(4),
        minor_axis_length=ones(4),
        major_axis_length=ones(4),
        perimeter=ones(4),
    )

    converted = converttounits(propdf, latlondata)

    top_left, top_right, bottom_left, bottom_right = 1, 2, 3, 4

    # For this georeferenced image, northern pixels (top rows) should have higher latitude.
    @test converted.latitude[top_left] > converted.latitude[bottom_left]
    @test converted.latitude[top_right] > converted.latitude[bottom_right]

    # Longitudes should increase from left to right at both top and bottom rows.
    @test converted.longitude[top_right] > converted.longitude[top_left]
    @test converted.longitude[bottom_right] > converted.longitude[bottom_left]
end

@testitem "converttounits copes with empty data frames" begin
    using IceFloeTracker
    import DataFrames: DataFrame

    refimage = "test_inputs/latlon/latlon_test_image-2020-06-21T00_00_00Z.tif"
    latlondata = latlon(refimage)
    latitude = latlondata[:latitude]
    nrows, ncols = size(latitude)

    propdf = DataFrame(;
        row_centroid=[],
        col_centroid=[],
        area=Float64[],
        convex_area=Float64[],
        minor_axis_length=Float64[],
        major_axis_length=Float64[],
        perimeter=Float64[],
    )

    @test converttounits(propdf, latlondata) isa DataFrame # doesn't crash
end

@testitem "_component_moment_measures: empty input" begin
    using IceFloeTracker.Segmentation: _component_moment_measures
    result = _component_moment_measures([0 0 0; 0 0 0; 0 0 0], [])
    @test result == Dict(
        :minor_axis_length=>[],
        :major_axis_length=>[],
        :row_centroid=>[],
        :col_centroid=>[],
        :orientation=>[],
    )
end

@testitem "_component_moment_measures: happy path" begin
    using IceFloeTracker.Segmentation: _component_moment_measures
    result = _component_moment_measures([0 0 0 0 0; 0 1 1 1 0; 0 0 0 0 0], [1])
    @test result == Dict(
        :minor_axis_length => [0.0],
        :major_axis_length => [3.265986323710904],
        :row_centroid => [2.0],
        :col_centroid => [3.0],
        :orientation => [1.5707963267948966],
    )
    result = _component_moment_measures([0 1 1 1 0; 0 2 2 2 0; 0 3 3 3 0], [1, 2, 3])
    @test result == Dict(
        :minor_axis_length => [0.0, 0.0, 0.0],
        :major_axis_length => [3.265986323710904, 3.265986323710904, 3.265986323710904],
        :row_centroid => [1.0, 2.0, 3.0],
        :col_centroid => [3.0, 3.0, 3.0],
        :orientation => [1.5707963267948966, 1.5707963267948966, 1.5707963267948966],
    )
end

@testitem "regionprops: bboxes: happy path" begin
    result = regionprops([0 0 0 0 0; 0 1 1 1 0; 0 0 0 0 0]; properties=[:bbox])
    @test result == Dict(:min_row => [2], :max_row => [2], :min_col => [2], :max_col => [4])

    result = regionprops([0 1 1 1 0; 0 2 2 2 0; 0 3 3 3 0]; properties=[:bbox])
    @test result == Dict(
        :min_col => [2, 2, 2],
        :max_col => [4, 4, 4],
        :min_row => [1, 2, 3],
        :max_row => [1, 2, 3],
    )
end

@testitem "bboxes for too-small regions" begin
    result = regionprops([0 0 0; 0 1 0; 0 0 0]; properties=[:bbox])
    @test result ==
        Dict(:min_row => Int[], :max_row => Int[], :min_col => Int[], :max_col => Int[])
end
@testitem "_count_pixels_in_hull: well-known shapes" begin
    import Images: convexhull
    import IceFloeTracker: PixelConvexArea
    import IceFloeTracker.Segmentation: _count_pixels_in_hull

    count_for(mask) = _count_pixels_in_hull(mask, convexhull(mask))

    @testset "filled rectangle: every pixel is inside its own hull" begin
        mask = trues(4, 6)
        @test count_for(mask) == 24
    end

    @testset "hollow rectangle: interior background lies inside the hull" begin
        mask = falses(5, 7)
        mask[1, :] .= true
        mask[end, :] .= true
        mask[:, 1] .= true
        mask[:, end] .= true
        # hull spans the full frame, so the hole is counted too
        @test count_for(mask) == 35
    end

    @testset "plus shape: hull is the diamond spanned by the arm tips" begin
        mask = falses(5, 5)
        mask[3, :] .= true
        mask[:, 3] .= true
        # pixels with |row - 3| + |col - 3| <= 2: 1 + 3 + 5 + 3 + 1
        @test count_for(mask) == 13
    end

    @testset "filled diamond: slanted edges exclude the corner background" begin
        mask = [abs(i - 3) + abs(j - 3) <= 2 for i in 1:5, j in 1:5]
        # bbox corners like (1,1) fail the diagonal edge tests, so only the
        # 13 diamond pixels are counted
        @test count_for(mask) == sum(mask) == 13
    end

    @testset "consistency with PixelConvexArea on a labeled image" begin
        labels = zeros(Int, 9, 9)
        labels[3:7, 3:7] .= 1  # 5x5 square
        @test PixelConvexArea()(labels)[1] == 25
    end

    @testset "KNOWN ISSUE: convexhull errors on right-triangle masks" begin
        # upstream ImageMorphology.convexhull throws "Not enough points to compute
        # convex hull" for staircase right triangles of any size, so their convex
        # area is reported as NaN via _convexhull_or_nothing
        mask = [i >= j for i in 1:6, j in 1:6]
        @test count_for(mask) == sum(mask) == 21 broken = true
        @test isnan(PixelConvexArea()(Int.(mask))[1])
    end
end

@testitem "convex-area guard marks background and sub-minimum labels undefined" begin
    import IceFloeTracker: PixelConvexArea, PolygonConvexArea

    # regression guard for the operator-precedence bug where the background label
    # fell through to a full convex-hull computation
    labels = zeros(Int, 9, 9)
    labels[3:7, 3:7] .= 1   # valid component, area 25: computed
    labels[9, 9] = 2        # too-small component, area 1 < minimum_area: undefined

    @testset "$(nameof(typeof(alg)))" for alg in (PixelConvexArea(), PolygonConvexArea())
        convex_areas = alg(labels)
        @test isnan(convex_areas[0]) # background: never a component
        @test isnan(convex_areas[2]) # sub-minimum component: undefined
        @test !isnan(convex_areas[1]) # valid component: computed
    end
end

@testitem "labeled-array validation rejects inputs without background-0 minimum" begin
    using IceFloeTracker.Segmentation: component_floes, component_convex_areas

    starts_at_two = [2 2; 2 2]
    has_negative = [0 1; -1 1]

    @test_throws ArgumentError component_floes(starts_at_two)
    @test_throws ArgumentError component_floes(has_negative)
    @test_throws ArgumentError component_convex_areas(starts_at_two)
    @test_throws ArgumentError component_convex_areas(has_negative)

    # minimum of 1 (no background pixels) is explicitly allowed
    @test component_floes([1 1; 1 1]) isa Dict
end

@testitem "PolygonConvexArea: single-argument functor on known shapes" begin
    using IceFloeTracker.Segmentation: PolygonConvexArea, PixelConvexArea

    @testset "square: shoelace over corner pixel centers" begin
        labels = zeros(Int, 9, 9)
        labels[3:7, 3:7] .= 1  # 5x5 square: hull corners span a 4x4 polygon
        @test PolygonConvexArea()(labels)[1] == 16.0
        # docstring invariant: polygon area is smaller than the pixel count
        @test PolygonConvexArea()(labels)[1] < PixelConvexArea()(labels)[1]
    end

    @testset "plus: diamond hull with both diagonals of length 4" begin
        plus = zeros(Int, 7, 7)
        plus[4, 2:6] .= 1
        plus[2:6, 4] .= 1
        @test PolygonConvexArea()(plus)[1] == 8.0  # d1 * d2 / 2
    end
end

@testitem "regionprops: gapped labelings report only the labels present" begin
    using IceFloeTracker: regionprops_table, PixelConvexArea
    import DataFrames: nrow
    import Images: component_lengths

    # A labeling with a gap in its numbering: labels 1 and 3 are present, 2 is
    # not. `label_components` never produces this, but dropping a region upstream
    # does, and `regionprops` takes the labeled array it is handed.
    labels = zeros(Int, 8, 8)
    labels[2:4, 2:4] .= 1
    labels[6:7, 6:7] .= 3

    # `component_lengths` counts every value in 0:maximum, so the absent label
    # gets a zero rather than no entry at all.
    areas = component_lengths(labels)
    @test axes(areas, 1) == 0:3
    @test areas[2] == 0

    # Every property family, at the default minimum_area and at 0. Single-pixel
    # regions survive minimum_area=0, an absent one still must not.
    for minimum_area in (1, 0)
        for properties in (
            [:label, :area],
            [:label, :perimeter],
            [:label, :circularity],
            [:label, :convex_area, :solidity],
            [:label, :bbox],
            [:label, :mask],
        )
            df = regionprops_table(
                labels;
                properties=properties,
                minimum_area=minimum_area,
                convex_area_algorithm=PixelConvexArea(),
            )
            @test df.label == [1, 3]
            @test nrow(df) == 2
        end
    end

    # Areas are the real pixel counts, not shifted by the gap.
    df = regionprops_table(labels; properties=[:label, :area])
    @test df.area == [9, 4]
end

@testsnippet IntensityFixture begin
    import StatsBase: mean, median, quantile
    import Images: RGB, SegmentedImage, segment_mean
    import IceFloeTracker: ExternalBoundary, expand_labels

    # Label 3 has area 1 and is dropped by the default `minimum_area=1`.
    labels = [
        1 1 0 0 0 0
        1 1 0 0 2 2
        0 0 0 0 2 2
        0 0 0 0 2 2
        3 0 0 0 0 0
    ]
    # Column-major pixel values 1:30, so label 1 holds 1, 2, 6, 7 and label 2 holds
    # 22, 23, 24, 27, 28, 29.
    intensity = Float64.(reshape(1:30, size(labels)))
    rgb_intensity = RGB.(intensity ./ 30, 0.5, 1 .- intensity ./ 30)

    # Helper function to call regionprops with the current labels and intensity arrays.
    regionprops(properties, extra_properties) = regionprops(
        labels, intensity; properties=properties, extra_properties=extra_properties
    )
end

@testitem "regionprops: extra_properties from plain functions" setup = [IntensityFixture] begin
    properties = [:label, :area]
    extra_properties = (mean, median)
    result = regionprops(properties, extra_properties)
    @test Set(keys(result)) ∈ Set([:label, :area, :mean, :median])
    @test result[:label] == [1, 2]
    @test result[:mean] == [4.0, 25.5]
    @test result[:median] == [4.0, 25.5]

    # vectors work as well as tuples
    props = regionprops([:label], [mean])
    @test props[:mean] == [4.0, 25.5]

    table = regionprops_table([:label], (maximum,)
    )
    @test names(table) ⊇ ["label", "maximum"]
    @test table.maximum == [7.0, 29.0]
end

@testitem "regionprops: extra_properties naming with pairs" setup = [IntensityFixture] begin
    extra_properties = (
        mean,
        :p90 => (v -> quantile(v, 0.9)),
        "smallest" => minimum,
        :wide_ring => ExternalBoundary(mean; radius=2),
    )
    result = regionprops([:label], extra_properties)
    @test Set(keys(result)) ⊆ Set([:label, :mean, :p90, :smallest, :wide_ring])
    @test result[:p90] ≈ [quantile([1.0, 2, 6, 7], 0.9), quantile([22.0, 23, 24, 27, 28, 29], 0.9)]
    @test result[:smallest] == [1.0, 22.0]

    # anonymous functions have no usable name
    @test_throws ArgumentError regionprops(
        labels, intensity; extra_properties=(v -> sum(v),)
    )
    # unsupported entries
    @test_throws ArgumentError regionprops(labels, intensity; extra_properties=(42,))
end

@testitem "regionprops: ExternalBoundary statistics" setup = [IntensityFixture] begin
    for r in (1, 2)
        result = regionprops(
            labels,
            intensity;
            properties=[:label],
            extra_properties=(ExternalBoundary(mean; radius=r),),
        )
        # Cross-check against the ring construction used in FSPipeline.
        ring = expand_labels(labels, r) .- labels
        expected = segment_mean(SegmentedImage(intensity, ring))
        @test result[:mean_boundary] ≈ [expected[1], expected[2]]
    end

    # the default radius is 15
    @test ExternalBoundary(mean).radius == 15
    @test_throws ArgumentError ExternalBoundary(mean; radius=0)

    # two radii side by side via pairs
    result = regionprops(
        labels,
        intensity;
        properties=[:label],
        extra_properties=(
            :ring1 => ExternalBoundary(mean; radius=1),
            :ring2 => ExternalBoundary(mean; radius=2),
        ),
    )
    @test result[:ring1] != result[:ring2]

    # no background around the floes: the ring is empty and the value is missing
    full = [1 1 2 2; 1 1 2 2]
    result = regionprops(
        full, rand(2, 4); properties=[:label], extra_properties=(ExternalBoundary(mean),)
    )
    @test all(ismissing, result[:mean_boundary])
end

@testitem "regionprops: extra_properties on color images" setup = [IntensityFixture] begin
    result = regionprops(labels, rgb_intensity; properties=[:label], extra_properties=(mean,))
    @test result[:mean][1] ≈ mean(rgb_intensity[labels .== 1])
    @test result[:mean][2] ≈ mean(rgb_intensity[labels .== 2])
    @test eltype(result[:mean]) <: RGB

    # A reducer that cannot handle the eltype (RGB has no ordering) raises its own error.
    @test_throws MethodError regionprops(labels, rgb_intensity; extra_properties=(median,))
    @test_throws MethodError regionprops(
        labels, rgb_intensity; extra_properties=(ExternalBoundary(median; radius=1),)
    )
end

@testitem "regionprops: extra_properties validation" setup = [IntensityFixture] begin
    # extras need an intensity image
    @test_throws ArgumentError regionprops(labels; extra_properties=(mean,))
    # sizes must match
    @test_throws DimensionMismatch regionprops(
        labels, intensity[1:(end-1), :]; extra_properties=(mean,)
    )
    # duplicated column names
    @test_throws ArgumentError regionprops(
        labels, intensity; extra_properties=(mean, :mean => median)
    )
    # collisions with requested base properties, including expanded ones
    @test_throws ArgumentError regionprops(
        labels, intensity; properties=[:area], extra_properties=(:area => mean,)
    )
    @test_throws ArgumentError regionprops(
        labels, intensity; properties=[:bbox], extra_properties=(:min_row => mean,)
    )
    # no extras: intensity image is still optional and unused
    @test regionprops(labels, intensity; properties=[:area]) ==
        regionprops(labels; properties=[:area])
end

@testitem "regionprops: extra_properties on an empty label image" setup = [IntensityFixture] begin
    import DataFrames: nrow

    empty_labels = zeros(Int, size(labels))
    table = regionprops_table(
        empty_labels,
        intensity;
        properties=[:label],
        extra_properties=(mean, ExternalBoundary(mean; radius=1)),
    )
    @test nrow(table) == 0
    @test names(table) ⊇ ["label", "mean", "mean_boundary"]
end
