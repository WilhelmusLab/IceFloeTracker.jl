#=
MWE: `PolygonConvexArea` produces `solidity > 1`, which is geometrically impossible.

Run (no dependencies beyond the package itself):

    julia --project=. benchmarks/issue_polygon_solidity_mwe.jl

Optionally append a labeled-scene fixture to add the real-data section:

    julia --project=. benchmarks/issue_polygon_solidity_mwe.jl --labels benchmarks/fixtures/labels_full.bin --rows 5680 --cols 3392

Solidity is `area / convex_area`. A region is always a subset of its own convex
hull, so solidity is bounded above by 1 and equals 1 exactly when the region is
already convex. With `convex_area_algorithm=PolygonConvexArea()` it exceeds 1
for every shape tested here, including a solid square.
=#

using IceFloeTracker: regionprops, PolygonConvexArea, PixelConvexArea

# ---------------------------------------------------------------- helpers

function argvalue(flag, default)
    i = findfirst(==(flag), ARGS)
    isnothing(i) && return default
    i == length(ARGS) && error("$flag requires a value")
    return ARGS[i + 1]
end

"""area, convex_area and solidity for a single-object label image."""
function measure(img, algorithm)
    d = regionprops(
        img; properties=[:area, :convex_area, :solidity], convex_area_algorithm=algorithm
    )
    return d[:area][1], d[:convex_area][1], d[:solidity][1]
end

function rule(title)
    println()
    println(title)
    println("-"^length(title))
end

# ---------------------------------------------------------------- 1. the MWE

rule("1. Minimal reproducer")

# A solid 3x3 square. It is convex, so it IS its own convex hull:
# convex_area must equal area (9) and solidity must be exactly 1.
square = zeros(Int64, 7, 7)
square[3:5, 3:5] .= 1

area, ca, sol = measure(square, PolygonConvexArea())
println("""
A solid 3x3 square (a convex region, so it is its own convex hull):

    square = zeros(Int64, 7, 7); square[3:5, 3:5] .= 1
    regionprops(square; properties=[:area, :convex_area, :solidity],
                convex_area_algorithm=PolygonConvexArea())

    area        = $area
    convex_area = $ca      <- expected $area
    solidity    = $sol   <- expected 1.0, and solidity > 1 is impossible
""")

# ---------------------------------------------------------------- 2. scope

rule("2. Not an edge case: every shape tested, including convex ones")

shapes = Pair{String,Matrix{Int64}}[
    "solid 3x3 square" => square,
    "solid 5x5 square" => (x = zeros(Int64, 9, 9); x[3:7, 3:7] .= 1; x),
    "solid 10x10 square" => (x = zeros(Int64, 14, 14); x[3:12, 3:12] .= 1; x),
    "disc r=5" => (
        x = zeros(Int64, 15, 15);
        for i in 1:15, j in 1:15
            (i - 8)^2 + (j - 8)^2 <= 25 && (x[i, j] = 1)
        end;
        x
    ),
    "L-shape" => (x = zeros(Int64, 11, 11); x[3:9, 3:4] .= 1; x[8:9, 3:9] .= 1; x),
    "T-shape" => (x = zeros(Int64, 11, 11); x[3:4, 3:9] .= 1; x[3:9, 6:7] .= 1; x),
]

println(rpad("shape", 20), rpad("area", 7), rpad("convex_area", 13),
        rpad("solidity", 10), "verdict")
for (name, img) in shapes
    a, c, s = measure(img, PolygonConvexArea())
    println(rpad(name, 20), rpad(a, 7), rpad(c, 13), rpad(round(s; digits=4), 10),
            s > 1 ? "IMPOSSIBLE" : "ok")
end

println("""

The first three are convex: their solidity must be exactly 1.0. The n x n square
gives n^2 / (n-1)^2, which tends to 1 only as n grows -- the error is worst on
the smallest regions, which is where floe segmentation produces the most objects.
""")

# ---------------------------------------------------------------- 3. cause

rule("3. Why")

println("""
`PolygonConvexArea` (src/Segmentation/regionprops.jl) applies the shoelace
formula to the convex hull returned by `ImageMorphology.convexhull`, whose
vertices are pixel *centres*:

    chull = _convexhull_or_nothing(A[bboxes[i]] .== i)
    ca = 0
    for j in 1:N
        x0, y0 = Tuple(chull[j])
        x1, y1 = Tuple(chull[(j%N)+1])
        ca += x0 * y1 - y0 * x1
    end
    ca *= 0.5

That is the area of the *polygon joining pixel centres*, a continuous area. But
`:area` is a pixel *count*. The two are not the same quantity, and the polygon
is always the smaller of the two.

Pick's theorem makes the gap exact. For a lattice polygon with `I` interior and
`B` boundary lattice points, its area is `A = I + B/2 - 1`, while the number of
pixels it covers is `I + B`. So

    pixel_count = polygon_area + B/2 + 1

The polygon area understates the covered pixel count by `B/2 + 1` -- roughly
half the hull perimeter -- for every region, with no exceptions. Dividing a
pixel count by it therefore always overshoots.
""")

a3, c3, _ = measure(square, PolygonConvexArea())
_, cp3, _ = measure(square, PixelConvexArea())
println("Checking that identity on the 3x3 square:")
println("  polygon_area                = ", c3)
println("  pixel count in hull         = ", cp3, "   (PixelConvexArea)")
println("  difference                  = ", cp3 - c3, "   = B/2 + 1 with B = 8")

# ---------------------------------------------------------------- 4. real data

rule("4. Scale on real data")

labelpath = argvalue("--labels", "")
if isempty(labelpath) || !isfile(labelpath)
    println("""
Skipped -- pass `--labels <file.bin> --rows R --cols C` to include this section.
Generate a fixture with `benchmarks/export_labels.jl`.
""")
else
    rows = parse(Int, argvalue("--rows", "0"))
    cols = parse(Int, argvalue("--cols", "0"))
    buf = Vector{Int32}(undef, rows * cols)
    read!(labelpath, buf)
    labels = Matrix{Int64}(permutedims(reshape(ltoh.(buf), cols, rows)))

    d = regionprops(
        labels; properties=[:area, :solidity], convex_area_algorithm=PolygonConvexArea()
    )
    s, a = d[:solidity], d[:area]
    finite = [i for i in eachindex(s) if isfinite(s[i])]
    bad = [i for i in finite if s[i] > 1]
    small = [i for i in finite if a[i] < 300]

    println("$(basename(labelpath)), $rows x $cols")
    println("  floes with finite solidity : ", length(finite))
    println("  solidity > 1 (impossible)  : ", length(bad),
            "  (", round(100 * length(bad) / length(finite); digits=1), "%)")
    println("  maximum solidity observed  : ", round(maximum(s[finite]); digits=4))
    println("  among floes < 300 px       : ", count(i -> s[i] > 1, small), " of ",
            length(small), " impossible")
end

# ---------------------------------------------------------------- 5. summary

rule("5. Summary")

println("""
Expected : solidity in (0, 1], equal to 1 for convex regions.
Actual   : solidity > 1 for every shape above, up to 2.25 on a solid 3x3 square.

Cause    : `:area` is a pixel count; `PolygonConvexArea` returns a continuous
           polygon area over pixel centres. `pixel_count = polygon_area + B/2 + 1`,
           so the denominator is systematically too small, worst for small regions.

`PixelConvexArea` (the default) counts pixels and does not have this problem --
it is only reachable via `convex_area_algorithm=PolygonConvexArea()`.

Secondary, not triggered in any case above: the shoelace sum has no `abs()`, so
its sign depends on the hull's vertex orientation. `ImageMorphology.convexhull`
appears to emit a consistent orientation -- no negative value was observed over
1819 real floes -- so this is a latent robustness issue rather than an active
bug, but it costs nothing to take the absolute value.
""")
