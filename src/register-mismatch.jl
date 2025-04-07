
"""
    mismatch(
        fixed::AbstractArray,
        moving::AbstractArray,
        test_angles::AbstractArray,
    )                   

Estimate a rotation that minimizes the 'mismatch' of aligning `moving` with `fixed`.

Returns a pair with the mismatch score `mm` and the associated registration angle `rot`.

# Arguments
- `fixed`,`moving`: images to align via a rigid transformation
- `test_angles`: candidate angles to check for rotations by, in degrees
```
"""
function mismatch(
    fixed::AbstractArray,
    moving::AbstractArray;
    test_angles=sort(reverse(range(; start=-180, stop=180, step=5)[1:(end-1)]); by=abs),
)
    shape_differences = shape_difference_rotation(fixed, moving, deg2rad.(test_angles))
    best_match = argmin((x) -> x.shape_difference, shape_differences)
    rotation_degrees = rad2deg(best_match.angle)
    normalized_area = (sum(fixed) + sum(moving)) / 2
    normalized_mismatch = best_match.shape_difference / normalized_area
    return (mm=normalized_mismatch, rot=rotation_degrees)
end
