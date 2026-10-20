import Images: imrotate, padarray, Fill
import Interpolations: BSpline, Constant
import StatsBase: mean

greaterthan05(x) = x .> 0.5 # used for the image resize step and for binarizing images
function imrotate_bin(x, r)
    return greaterthan05(collect(imrotate(x, r, axes(x); method=BSpline(Constant()))))
end
function imrotate_bin_nocrop(x, r)
    return greaterthan05(collect(imrotate(x, r; method=BSpline(Constant()))))
end
imrotate_bin_clockwise_radians(x, r) = imrotate_bin(x, r)
imrotate_bin_counterclockwise_radians(x, r) = imrotate_bin(x, -r)
imrotate_bin_clockwise_degrees(x, r) = imrotate_bin_clockwise_radians(x, deg2rad(r))
function imrotate_bin_counterclockwise_degrees(x, r)
    return imrotate_bin_counterclockwise_radians(x, deg2rad(r))
end

"""
    compute_centroid(im::AbstractArray{Bool}; rounded=false)

Calculate the centroid of a binary image. If 'rounded', return the
nearest integer.
"""
function compute_centroid(im::AbstractArray{Bool}; rounded=false)
    xi = 0
    yi = 0
    R = count(im)
    for idx in CartesianIndices(im)
        if im[idx]
            ii, jj = Tuple(idx)
            xi += ii
            yi += jj
        end
    end

    x0, y0 = xi / R, yi / R
    rounded && return round(Int32, x0), round(Int32, y0)
    return x0, y0
end

# floor-rounded centroid, i.e. the pixel containing the centroid
floor_centroid(im) = Int64.(floor.(compute_centroid(im; rounded=false)))

# symmetric difference of two equal-size binary masks, without temporaries
count_symdiff(im1, im2) = count(((a, b),) -> a ⊻ b, zip(im1, im2))

"""
    align_centroids(im1::AbstractArray{Bool}, im2::AbstractArray{Bool})
    align_centroids(im1, centroid1, im2, centroid2)

Align images by padding so that the centroids of each image are on the edge of or within the same pixel.
The 4-argument form takes precomputed [`floor_centroid`](@ref)s, so callers can reuse a
centroid across many alignments against the same image.
"""
function align_centroids(im1::AbstractArray{Bool}, im2::AbstractArray{Bool})
    return align_centroids(im1, floor_centroid(im1), im2, floor_centroid(im2))
end

function align_centroids(
    im1::AbstractArray{Bool},
    centroid1::Tuple{Int64,Int64},
    im2::AbstractArray{Bool},
    centroid2::Tuple{Int64,Int64},
)
    # Location of the pixel containing the centroids of im1 and im2
    # in their current coordinate systems
    r1, c1 = centroid1
    r2, c2 = centroid2

    # Calculate the same centroid, but measured from the bottom right of each image
    s1, d1 = size(im1) .- (r1, c1) .+ 1
    s2, d2 = size(im2) .- (r2, c2) .+ 1

    # Calculate the new "common centroid" position in image coordinates
    rn, cn = (max(r1, r2), max(c1, c2))
    # Calculate the new "reverse common centroid" position in image coordinates from the bottom right
    sn, dn = (max(s1, s2), max(d1, d2))

    # For each image, we shift the pixel containing its centroid to the new centroid
    # by adding rn-ri rows padding at the top, and cn-ci columns at the left.
    # We ensure that the centroid is the same distance from the right border
    # by adding sn-si rows padding at the bottom and dn-di columns padding at the right
    # These need to be `collect`
    im1_padded = collect(padarray(im1, Fill(0, (rn - r1, cn - c1), (sn - s1, dn - d1))))
    im2_padded = collect(padarray(im2, Fill(0, (rn - r2, cn - c2), (sn - s2, dn - d2))))

    return im1_padded, im2_padded
end

"""
    shape_difference(floe1::DataFrameRow, floe2::DataFrameRow)
    shape_difference(floe1_mask::BitMatrix, floe1_orientation::Float64,
                     floe2_mask::BitMatrix, floe2_orientation::Float64)

    Computes the shape difference between two ice floes using the estimated floe
    orientation to first align both floes on the same axis.

"""
function shape_difference(
    floe1_mask::AbstractArray{Bool},
    floe1_orientation::Real,
    floe2_mask::AbstractArray{Bool},
    floe2_orientation::Real,
)
    im1 = imrotate_bin_nocrop(floe1_mask, floe1_orientation)
    im2 = imrotate_bin_nocrop(floe2_mask, floe2_orientation)
    im1, im2 = align_centroids(im1, im2)
    return count_symdiff(im1, im2)
end

function shape_difference(floe1::DataFrameRow, floe2::DataFrameRow)
    return shape_difference(floe1.mask, floe1.orientation, floe2.mask, floe2.orientation)
end

"""
     shape_difference_rotation(im_reference, im_target, test_angles; imrotate_function=imrotate_bin_clockwise_radians)

Computes the shape difference between im_reference and im_target for each angle in test_angles.
The reference image is held constant, while the target image is rotated. The test_angles are interpreted
as the angle of rotation from target to reference, so to find the best match, we rotate the reverse
direction. A perfect match at angle `A` would imply im_target is the same shape as if im_reference was
rotated by `A`.
Use `imrotate_function=imrotate_bin_<clockwise|counterclockwise>_<radians|degrees>` to get angles <clockwise|counterclockwise> in <radians|degrees>.
"""
function shape_difference_rotation(
    im_reference, im_target, test_angles; imrotate_function=imrotate_bin_clockwise_radians
)
    shape_differences = Array{
        NamedTuple{(:angle, :shape_difference),Tuple{Float64,Float64}}
    }(
        undef, length(test_angles)
    )

    # the reference image never changes, so compute its centroid once
    centroid_reference = floor_centroid(im_reference)

    for (idx, angle) in enumerate(test_angles)

        # rotate image back by angle
        im_target_rotated = imrotate_function(im_target, -angle)

        im1, im2 = align_centroids(
            im_reference,
            centroid_reference,
            im_target_rotated,
            floor_centroid(im_target_rotated),
        )

        shape_difference = count_symdiff(im1, im2)
        shape_differences[idx] = (; angle, shape_difference)
    end
    return shape_differences
end

"""
The default registration angles are evenly distributed in steps of π/36 rad (5º) around a full rotation,
ensuring that no angles are repeated (since -π rad == π rad),
and ordered so that smaller absolute angles which are positive will be returned in the event of a tie in the shape difference.
"""
register_default_angles_rad = sort(
    reverse(range(; start=(-π), stop=π, step=π / 36)[1:(end - 1)]); by=abs
)
"""
Coarse registration angles for [`register_coarse_to_fine`](@ref): -π/6 to π/6 rad (±30º) in steps
of π/36 rad (5º), ordered like `register_default_angles_rad` so that smaller absolute angles which are
positive are returned in the event of a tie in the shape difference.
"""
register_coarse_angles_rad = sort(
    range(; start=(-π / 6), stop=π / 6, step=π / 36); by=θ -> (abs(θ), -θ)
)

"""
    register_coarse_to_fine(
        im_reference,
        im_target;
        coarse_angles=register_coarse_angles_rad,
        fine_step=deg2rad(1.0),
        fine_window=deg2rad(4.0),
        imrotate_function=imrotate_bin_clockwise_radians,
    )

Finds the rotation angle between `im_reference` and `im_target` in two stages: [`register`](@ref) over
`coarse_angles`, then over the angles within `fine_window` of the coarse result in steps of `fine_step`.
`fine_step` must be finite and positive, and `fine_window` finite and non-negative; a `fine_window` of 0
gives a coarse-only search. Throws an `ArgumentError` otherwise.
Ties are broken toward smaller absolute angles which are positive, as in `register`.
Logs a warning when the coarse result is the smallest or largest angle in `coarse_angles`, since the
true rotation may then lie outside the searched range; `coarse_angles` is treated as an interval, so a
full-circle grid warns at its ends.
Angles follow the convention of `imrotate_function`, in radians.
"""
function register_coarse_to_fine(
    im_reference,
    im_target;
    coarse_angles=register_coarse_angles_rad,
    fine_step::Real=deg2rad(1.0),
    fine_window::Real=deg2rad(4.0),
    imrotate_function=imrotate_bin_clockwise_radians,
)
    0 < fine_step < Inf ||
        throw(ArgumentError("fine_step must be finite and positive, got $fine_step"))
    0 <= fine_window < Inf || throw(
        ArgumentError("fine_window must be finite and non-negative, got $fine_window")
    )
    coarse = register(im_reference, im_target; test_angles=coarse_angles, imrotate_function)
    lo, hi = extrema(coarse_angles)
    if coarse ≈ lo || coarse ≈ hi
        @warn "coarse registration angle $(round(rad2deg(coarse); digits=1))° is on the edge of the search range"
    end
    n = floor(Int, fine_window / fine_step)
    fine_angles = sort([coarse + k * fine_step for k in (-n):n]; by=θ -> (abs(θ), -θ))
    return register(im_reference, im_target; test_angles=fine_angles, imrotate_function)
end

"""
    register(
        im_reference,
        im_target;
        test_angles=register_default_angles_rad,
        imrotate_function=imrotate_bin_clockwise_radians,
    )

Finds the image rotation angle in `test_angles` which minimizes the shape difference between `im_reference` and `im_target`.
The default test angles are shown in `register_default_angles_rad`.
Use `imrotate_function=imrotate_bin_<clockwise|counterclockwise>_<radians|degrees>` to get angles <clockwise|counterclockwise> in <radians|degrees>.
"""
function register(
    im_reference,
    im_target;
    test_angles=register_default_angles_rad,
    imrotate_function=imrotate_bin_clockwise_radians,
)
    shape_differences = shape_difference_rotation(
        im_reference, im_target, test_angles; imrotate_function
    )
    best_match = argmin((x) -> x.shape_difference, shape_differences)
    return best_match.angle
end

"""
    mismatch(
        fixed::AbstractArray,
        moving::AbstractArray,
        test_angles::AbstractArray,
    )

Estimate a rotation that minimizes the 'mismatch' of aligning `moving` with `fixed`.

Returns a pair with the mismatch score `mm` and the associated registration angle `rot`.

## Arguments
- `fixed`,`moving`: images to align via a rigid transformation
- `test_angles`: candidate angles to check for rotations by, in degrees.
  In the case of a tie in the shape difference, the earlier angle from this array will be returned.
"""
function mismatch(fixed::AbstractArray, moving::AbstractArray, test_angles::AbstractArray)
    shape_differences = shape_difference_rotation(
        fixed, moving, test_angles; imrotate_function=imrotate_bin_clockwise_degrees
    )
    best_match = argmin((x) -> x.shape_difference, shape_differences)
    rotation_degrees = best_match.angle
    normalized_area = (sum(fixed) + sum(moving)) / 2
    normalized_mismatch = best_match.shape_difference / normalized_area
    return (mm=normalized_mismatch, rot=rotation_degrees)
end

"""
    mismatch(
        fixed::AbstractArray,
        moving::AbstractArray,
        mxrot::Real,
        step::Real,
    )

Estimate a rotation that minimizes the 'mismatch' of aligning `moving` with `fixed`.

Returns a pair with the mismatch score `mm` and the associated registration angle `rot`.

## Arguments
- `fixed`,`moving`: images to align via a rigid transformation
- `mxrot`: maximum rotation angle in degrees
- `step`: rotation angle step size in degrees

The default registration angles are evenly distributed in steps of 5º around a full rotation,
ensuring that no angles are repeated (since -180º == +180º).

Angles are ordered so that smaller absolute angles which are positive will be returned in the event of a tie in the shape difference.
"""
function mismatch(
    fixed::AbstractArray, moving::AbstractArray, mxrot::Real=180, step::Real=5
)
    test_angles = sort(
        reverse(range(; start=(-mxrot), stop=mxrot, step=step)[1:(end - 1)]); by=abs
    )
    return mismatch(fixed, moving, test_angles)
end

# ============================================================================
# Geometric Transformations for Boundary Curves
# ============================================================================

"""
    _get_rotation_matrix(angle::Real)

Create a 2D rotation matrix for the given angle in radians.
Positive angle = counterclockwise rotation.
"""
function _get_rotation_matrix(angle::Real)
    cos_a = cos(angle)
    sin_a = sin(angle)
    return [cos_a -sin_a; sin_a cos_a]
end

# The points of a curve without the closing duplicate. Traced boundaries repeat their
# first point as the last (after resampling, only up to round-off of order 1e-13), so the
# duplicate would be counted twice in a centroid and make it depend on the trace start.
# Open curves are returned unchanged.
function _normalize_contour(boundary::AbstractMatrix)
    n = size(boundary, 1)
    is_closed = n > 1 && isapprox(@view(boundary[1, :]), @view(boundary[n, :]); atol=1e-9)
    return is_closed ? @view(boundary[1:(end-1), :]) : boundary
end

"""
    rotate_boundary(boundary::Matrix{Float64}, angle::Real; center::Union{Nothing,Tuple{Float64,Float64}}=nothing)

Apply 2D rotation to boundary curve around center point using matrix multiplication.
Angle is in radians, positive = counterclockwise.

# Arguments
- `boundary`: Matrix(n, 2) with [x y] coordinates
- `angle`: Rotation angle in radians
- `center`: Center of rotation; if nothing, uses centroid of boundary
"""
function rotate_boundary(
    boundary::Matrix{Float64},
    angle::Real;
    center::Union{Nothing,Tuple{Float64,Float64}}=nothing,
)
    center = isnothing(center) ? vec(mean(_normalize_contour(boundary); dims=1)) : collect(center)
    rot_matrix = _get_rotation_matrix(angle)
    boundary_centered = boundary .- center'
    rotated_centered = boundary_centered * transpose(rot_matrix)
    return rotated_centered .+ center'
end

"""
    center_boundary(boundary::Matrix{Float64}; target_center::Tuple{Float64,Float64}=(0.0, 0.0))

Translate boundary curve to center at target_center.

# Arguments
- `boundary`: Matrix(n, 2) with [x y] coordinates
- `target_center`: Target centroid position (default: origin)
"""
function center_boundary(
    boundary::Matrix{Float64}; target_center::Tuple{Float64,Float64}=(0.0, 0.0)
)
    centroid = vec(mean(_normalize_contour(boundary); dims=1))
    offset = collect(target_center) .- centroid
    return boundary .+ offset'
end

# ============================================================================
# Distance Metrics for Boundary Curves
# ============================================================================
#
# Both metrics treat a boundary as a point SET, so they do not depend on where the trace
# happens to start. An index-wise comparison (point i against point i) does: two traces
# of the same floe begin at unrelated points, and on such input an index-wise metric
# recovered the true rotation angle on 2% of real floes.

# Nearest-neighbour distance from point `a` to the point set `B`. Compared on squared
# distances; the square root is taken once.
function _nn_dist(a, B)
    best = Inf
    @inbounds for i in axes(B, 1)
        d2 = (a[1] - B[i, 1])^2 + (a[2] - B[i, 2])^2
        d2 < best && (best = d2)
    end
    return sqrt(best)
end

"""
    boundary_hausdorff(b1::Matrix{Float64}, b2::Matrix{Float64})

Hausdorff distance between two boundaries treated as point sets, after centering both at
the origin: the larger of the two directed distances, each the greatest nearest-neighbour
distance from one set to the other. Insensitive to where either trace begins. Units:
pixels.

Reference: Huttenlocher, Klanderman & Rucklidge (1993), *Comparing images using the
Hausdorff distance*, IEEE TPAMI 15(9).
"""
function boundary_hausdorff(b1::Matrix{Float64}, b2::Matrix{Float64})
    A = center_boundary(Matrix(_normalize_contour(b1)))
    B = center_boundary(Matrix(_normalize_contour(b2)))
    directed(P, Q) = maximum(_nn_dist(@view(P[i, :]), Q) for i in axes(P, 1))
    return max(directed(A, B), directed(B, A))
end

"""
    boundary_modified_hausdorff(b1::Matrix{Float64}, b2::Matrix{Float64})

Modified Hausdorff distance between two boundaries treated as point sets, after centring
both at the origin. Each directed distance is the *mean* rather than the maximum
nearest-neighbour distance, so a single outlying point does not set the score.
Insensitive to where either trace begins. Units: pixels.

Reference: Dubuisson & Jain (1994), *A modified Hausdorff distance for object matching*,
Proc. 12th IAPR Int. Conf. on Pattern Recognition.
"""
function boundary_modified_hausdorff(b1::Matrix{Float64}, b2::Matrix{Float64})
    A = center_boundary(Matrix(_normalize_contour(b1)))
    B = center_boundary(Matrix(_normalize_contour(b2)))
    directed(P, Q) = sum(_nn_dist(@view(P[i, :]), Q) for i in axes(P, 1)) / size(P, 1)
    return max(directed(A, B), directed(B, A))
end
