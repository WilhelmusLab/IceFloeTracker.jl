# `PolygonConvexArea` returns `solidity > 1`

## Summary

`regionprops(...; convex_area_algorithm=PolygonConvexArea())` returns
`:solidity` values greater than 1. Solidity is `area / convex_area`, and a
region is always a subset of its own convex hull, so it is bounded above by 1
and equals 1 exactly when the region is already convex.

It happens for **every** shape tested, including convex ones, and for **57.4%
of the 1819 floes** in a real Fram Strait scene (max observed 3.33).

The default algorithm, `PixelConvexArea`, is unaffected.

**This is not a newly discovered approximation — see "What is already
documented" below.** The polygon method's undercount is described in the
existing docstrings. What appears not to have been considered is that dividing
a pixel count by it, which `:solidity` now does, yields values outside the
mathematically possible range.

## Reproducer

```
julia --project=. benchmarks/issue_polygon_solidity_mwe.jl
```

A solid 3×3 square is convex, so it *is* its own convex hull:

```julia
using IceFloeTracker: regionprops, PolygonConvexArea

square = zeros(Int64, 7, 7); square[3:5, 3:5] .= 1
regionprops(square; properties=[:area, :convex_area, :solidity],
            convex_area_algorithm=PolygonConvexArea())
```

```
area        = 9
convex_area = 4.0     # expected 9
solidity    = 2.25    # expected 1.0
```

| shape | area | convex_area | solidity | |
|---|---:|---:|---:|---|
| solid 3x3 square | 9 | 4.0 | 2.25 | impossible |
| solid 5x5 square | 25 | 16.0 | 1.5625 | impossible |
| solid 10x10 square | 100 | 81.0 | 1.2346 | impossible |
| disc r=5 | 81 | 74.0 | 1.0946 | impossible |
| L-shape | 24 | 23.5 | 1.0213 | impossible |
| T-shape | 24 | 23.5 | 1.0213 | impossible |

The first three are convex and must give exactly 1.0. An n×n square gives
`n²/(n-1)²`, which approaches 1 only as n grows — so the error is worst on the
smallest regions, which is where floe segmentation produces the most objects.

## What is already documented

The direction and size-dependence of the discrepancy are already stated in the
code. `component_convex_areas`:

> The polygon method uses Green's theorem to find the area of a polygon through
> its line integral, while the pixel method uses a point-in-pixel calculation to
> determine if pixels are inside the convex hull. **In general the polygon area
> will be smaller than the pixel area.**

and `PolygonConvexArea`:

> Estimate the convex area by integrating the area of the convex hull polygon.
> [...] **In general, the error should be smaller for larger shapes.**

The Green's theorem approach was the deliberate original design, proposed in
#820. So the approximation itself is known and intended.

What changed is that `:solidity` was added later (`df398d4c`, 2026-09-11),
after the convex-area algorithms (`c25476ef`, 2025-12-19), and it divides
`:area` — a pixel count — by `:convex_area`. Under `PolygonConvexArea` those
are two different kinds of quantity, so the documented "polygon area is
smaller" becomes a ratio above 1. A search of the issue tracker for
"solidity" and "convex" turned up no existing report of this.

The question this issue raises is therefore not "is the polygon area
approximate" — that is documented — but "should a property with a hard
mathematical bound of 1 be allowed to return 2.25".

## Cause

`PolygonConvexArea` (`src/Segmentation/regionprops.jl`) applies the shoelace
formula to the hull from `ImageMorphology.convexhull`, whose vertices are pixel
**centres**:

```julia
chull = _convexhull_or_nothing(A[bboxes[i]] .== i)
ca = 0
for j in 1:N
    x0, y0 = Tuple(chull[j])
    x1, y1 = Tuple(chull[(j%N)+1])
    ca += x0 * y1 - y0 * x1
end
ca *= 0.5
```

That is the area of the polygon joining pixel centres — a continuous area.
But `:area` is a pixel **count**. They are different quantities, and the
polygon is always the smaller one.

Pick's theorem makes the gap exact. For a lattice polygon with `I` interior and
`B` boundary lattice points, the area is `A = I + B/2 - 1` while the pixels it
covers number `I + B`, so

```
pixel_count = polygon_area + B/2 + 1
```

The denominator is short by `B/2 + 1` — about half the hull perimeter — for
every region without exception, so the ratio always overshoots. On the 3×3
square: polygon area 4.0, pixels in hull 9, difference 5 = 8/2 + 1.

## Scale on real data

`001-fram_strait-20120412.aqua.labeled.png`, 5680×3392:

```
floes with finite solidity : 1819
solidity > 1 (impossible)  : 1045  (57.4%)
maximum solidity observed  : 3.3333
among floes < 300 px       : 862 of 1040 impossible
```

## Suggested fix

Options, in rough order of increasing scope:

1. **Make `:solidity` refuse the mismatch.** The cheapest correct move is to
   compute `:solidity` from a pixel-counted convex area regardless of
   `convex_area_algorithm`, so the ratio always compares like with like.
2. **Document the restriction** — state in `PolygonConvexArea`'s docstring that
   it is not suitable as a `:solidity` denominator, per the repo's convention
   of putting known defects in a comment at the code concerned.
3. **Fix the underlying quantity**, below, which also resolves a separate
   disagreement with scikit-image.

Make the hull enclose the pixel *extents* rather than their centres, which is
what scikit-image's `convex_hull_image` does by default (`offset_coordinates=True`:
each pixel contributes four points at ±0.5 along each axis). Hulling those and
then counting covered pixels makes `convex_area` commensurable with `area`
again, and makes `PixelConvexArea` agree with `skimage.measure`'s `area_convex`
— which it currently does not, for the same underlying reason. (Verified
separately: Julia matches `convex_hull_image(..., offset_coordinates=False)`
exactly on 27/27 floes, and the default on 0/27.)

The offsets can be handled in exact integer arithmetic by doubling coordinates:
pixel `(r,c)` contributes `(2r±1, 2c)` and `(2r, 2c±1)`, and pixel centres test
as `(2r, 2c)`. That needs a convex hull over arbitrary integer points;
`ImageMorphology.convexhull` only accepts a binary image, so it cannot be
reused directly.

Secondary, latent rather than active: the shoelace sum has no `abs()`, so its
sign depends on the hull's vertex orientation. No negative value appeared over
the 1819 real floes, so `convexhull` appears to emit a consistent orientation,
but taking the absolute value costs nothing.
