# `convex_area` disagrees with scikit-image: cause confirmed

Found while building `benchmark_regionprops.jl`. Recorded here rather than in a
docstring, per the repo's docstring conventions.

## What happens

On the 512x512 Fram Strait crop, `IceFloeTracker.regionprops`'s `:convex_area`
is **lower than** `skimage`'s `area_convex` for **27 of 27** floes -- never
higher, never equal. Relative gap 1.0% to 12.5%, median 6.2%, rising to 31% on
the full scene, and worse for smaller floes. `:solidity` inherits all of it.

## Cause: the hull encloses pixel centres, not pixel extents

`skimage.morphology.convex_hull_image` defaults to `offset_coordinates=True`,
which replaces each pixel `(r, c)` with four points at `(r±0.5, c)` and
`(r, c±0.5)` before hulling. The hull is therefore inflated by half a pixel
outward. IceFloeTracker hulls the bare pixel centres:
`ImageMorphology.convexhull` takes a binary image and returns
`Vector{CartesianIndex{2}}`, integer coordinates only.

That single difference accounts for all of it:

| | T-shape | 512 crop, 27 floes |
|---|---|---|
| IceFloeTracker | 30 | -- |
| skimage, `offset_coordinates=False` | 30 | **27/27 exact match with Julia** |
| skimage, default (`offset_coordinates=True`) | 34 | 0/27 match |

Reproduce the second column with:

```python
convex_hull_image(p.image, offset_coordinates=False).sum()
```

## A hypothesis that was wrong

An earlier version of this file blamed a strict inequality excluding pixels
lying exactly on a hull edge. That is not it. `_count_pixels_in_hull` rejects
on `edge_cross_product < 0`, so on-edge pixels *are* counted -- which already
matches skimage's `include_borders=True` (`labels >= 1`, counting vertex and
edge classifications as inside). Border handling agrees; only the half-pixel
offset differs.

## Minimal reproducer

```julia
using IceFloeTracker: regionprops
T = zeros(Int64, 11, 11); T[3:4, 3:9] .= 1; T[3:9, 6:7] .= 1
regionprops(T; properties=[:area, :convex_area])
# area = 24, convex_area = 30.0
```

```python
import numpy as np
from skimage.measure import regionprops_table
T = np.zeros((11, 11), dtype=np.int32); T[2:4, 2:9] = True; T[2:9, 5:7] = True
regionprops_table(T, properties=('area', 'area_convex'))
# area = 24, area_convex = 34.0
```

Convex shapes agree trivially (a solid 5x5 square gives 25 both sides), because
a half-pixel dilation of a shape that already fills its hull adds no new pixel
centres. The gap appears only where the hull cuts diagonally across background.

## What matching skimage would require

Hull the offset points, then count pixel centres inside. This can be done in
exact integer arithmetic by doubling coordinates -- pixel `(r,c)` contributes
`(2r±1, 2c)` and `(2r, 2c±1)`, centres test as `(2r, 2c)` -- avoiding
floating-point robustness questions. It needs a convex hull over arbitrary
integer points; `ImageMorphology.convexhull` only accepts a binary image, so it
cannot be reused directly.

Two consequences worth noting before doing it:

- `:convex_area` feeds `convex_area_relative_error_filter`, which is in the
  **default** tracking filter set (`src/Tracking/filter_functions.jl:236,272`,
  mirrored in `src/Pipeline/FSPipeline.jl:723,759`). It is a relative-error
  filter, so a systematic inflation on both sides of a pair largely cancels,
  but the thresholds were calibrated against current values.
- Offset hulls are well defined for regions of any size, so the
  `minimum_area = 4` guard and the hull-failure `NaN` become unnecessary
  (skimage returns `convex_area = 1, solidity = 1` for a single pixel). That
  would remove the 15 non-comparable rows on the full scene, and would change
  the existing `@test isnan(...)` for the 4-pixel cross in
  `test/test-regionprops.jl:102`.

## Separate defect in the other algorithm

`PolygonConvexArea` has a distinct and more serious problem -- it returns
`solidity > 1` for 57% of real floes. See `issue_polygon_solidity.md` and
`issue_polygon_solidity_mwe.jl`.
