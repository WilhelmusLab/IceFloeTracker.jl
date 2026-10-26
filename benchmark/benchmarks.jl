using IceFloeTracker, BenchmarkTools
using DataFrames: DataFrame
using FileIO: load
using Images: Gray, label_components

const SUITE = BenchmarkGroup()

# Inputs come from tracked test data via `pkgdir`, so each benchmarked revision reads
# its own checkout regardless of where the script runs from.
const SCENE = pkgdir(IceFloeTracker, "test", "test_inputs", "matlab_isolated_floes.png")
const LABELS = Matrix{Int64}(label_components(Gray.(load(SCENE)) .> 0))

# =============================================================================
# regionprops
# =============================================================================

SUITE["regionprops"] = BenchmarkGroup()

for prop in (:area, :perimeter, :convex_area, :bbox, :centroid, :orientation)
    SUITE["regionprops"][string(prop)] = @benchmarkable regionprops(
        $LABELS; properties=[$prop], convex_area_algorithm=PixelConvexArea()
    )
end
SUITE["regionprops"]["default_properties"] = @benchmarkable regionprops(
    $LABELS; convex_area_algorithm=PixelConvexArea()
)

# =============================================================================
# Floe masks and ψ-s curves
# =============================================================================

SUITE["extend_regionprops"] = BenchmarkGroup()

const PROPS = DataFrame(regionprops(LABELS; properties=[:label, :area]))

SUITE["extend_regionprops"]["add_floemasks!"] = @benchmarkable add_floemasks!(
    df, $LABELS
) setup = (df = copy($PROPS))

const PROPS_MASKED = let df = copy(PROPS)
    add_floemasks!(df, LABELS)
    df
end

SUITE["extend_regionprops"]["add_ψs!"] = @benchmarkable add_ψs!(df) setup = (
    df = copy($PROPS_MASKED)
)

# =============================================================================
# Registration
# =============================================================================

SUITE["register"] = BenchmarkGroup()

# Pad, rotate, then re-crop to the bounding box, so no part of the floe is clipped.
function rotated_mask(m, θ)
    p = ceil(Int, 0.5 * maximum(size(m)))
    padded = falses(size(m) .+ 2p)
    padded[(p + 1):(p + size(m, 1)), (p + 1):(p + size(m, 2))] .= m
    r = imrotate_bin_counterclockwise_radians(padded, θ)
    rows = findall(vec(any(r; dims=2)))
    cols = findall(vec(any(r; dims=1)))
    return r[first(rows):last(rows), first(cols):last(cols)]
end

# One floe per area bucket: the median-area floe inside the bucket.
for (lo, hi) in ((100, 300), (300, 1000), (1000, 5000), (5000, typemax(Int)))
    in_bucket = findall(a -> lo <= a < hi, PROPS_MASKED.area)
    isempty(in_bucket) && continue
    i = in_bucket[sortperm(PROPS_MASKED.area[in_bucket])[cld(length(in_bucket), 2)]]
    mask = PROPS_MASKED.mask[i]
    target = rotated_mask(mask, deg2rad(30))
    name = hi == typemax(Int) ? "area_$(lo)+" : "area_$(lo)-$(hi)"
    SUITE["register"][name] = @benchmarkable register($mask, $target)
end
