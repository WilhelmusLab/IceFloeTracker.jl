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

    # Integer padding shifts both floor-centroids exactly to (rn, cn), so they stay aligned
    @assert size(im1_padded) == size(im2_padded)

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

        # Check here that im1 and im2 sizes are the same
        # This should be guaranteed by "align_centroids"
        @assert size(im1) == size(im2)

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
# normalize to [-π, π), the convention of register_default_angles_rad
function normalize_angle(θ)
    θn = rem2pi(θ, RoundNearest)
    return θn == π ? oftype(θn, -π) : θn
end

"""
    prior_test_angles(prior_rad; window=deg2rad(10.0), step=π / 180)

Build registration test angles concentrated around a prior rotation estimate `prior_rad`
(in radians, in the convention of `imrotate_bin_clockwise_radians`) and its 180° alias
`prior_rad + π`, since priors derived from image-moment orientations are only defined modulo π.
Offsets are built symmetrically from zero to ensure the prior angle itself is always included,
even when the window is not an exact multiple of the step. Angles are normalized to [-π, π)
and ordered so that, in the event of a tie in the shape difference, smaller absolute angles
which are positive are preferred, matching `register_default_angles_rad`.
"""
function prior_test_angles(prior_rad::Real; window::Real=deg2rad(10.0), step::Real=π / 180)
    max_steps = ceil(Int, window / step)
    offsets = collect((-max_steps):max_steps) .* step
    offsets = filter(o -> abs(o) <= window, offsets)
    angles = [
        normalize_angle(alias + offset) for alias in (prior_rad, prior_rad + π) for
        offset in offsets
    ]
    unique!(angles)
    return sort!(angles; by=x -> (abs(x), -x))
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
    center = isnothing(center) ? vec(mean(boundary; dims=1)) : collect(center)
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
    centroid = vec(mean(boundary; dims=1))
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

# The point set of a boundary. Boundaries are closed (first row repeated as last), and
# which vertex carries the duplicate depends on where the trace started, so a mean over
# all rows is start-dependent; over the unique vertices it is not.
function _unique_points(b::Matrix{Float64})
    n = size(b, 1)
    return n > 1 && @views(b[1, :] == b[n, :]) ? b[1:(n - 1), :] : b
end

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

Hausdorff distance between two boundaries treated as point sets, after centring both at
the origin: the larger of the two directed distances, each the greatest nearest-neighbour
distance from one set to the other. Insensitive to where either trace begins. Units:
pixels.

Reference: Huttenlocher, Klanderman & Rucklidge (1993), *Comparing images using the
Hausdorff distance*, IEEE TPAMI 15(9).
"""
function boundary_hausdorff(b1::Matrix{Float64}, b2::Matrix{Float64})
    A = center_boundary(_unique_points(b1))
    B = center_boundary(_unique_points(b2))
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
    A = center_boundary(_unique_points(b1))
    B = center_boundary(_unique_points(b2))
    directed(P, Q) = sum(_nn_dist(@view(P[i, :]), Q) for i in axes(P, 1)) / size(P, 1)
    return max(directed(A, B), directed(B, A))
end

# ============================================================================
# Boundary-Curve Registration
# ============================================================================

"""
    shape_difference_rotation_boundary(boundary_reference, boundary_target, test_angles;
                                      metric=boundary_modified_hausdorff)

Boundary-curve analogue of [`shape_difference_rotation`](@ref). Computes the shape
difference between `boundary_reference` and `boundary_target` for each angle in
`test_angles`, holding the reference fixed and rotating the target.

Angle convention matches the mask-based version: a returned angle `A` means the target
mask looks like the reference mask rotated by `A` under `imrotate_bin_clockwise_radians`,
so `register_boundary` on traced boundaries agrees in sign with `register` on the masks
they were traced from.

`metric(reference, rotated_target)` may be any function returning a real shape
difference; see `boundary_modified_hausdorff` (default) and `boundary_hausdorff`.
"""
function shape_difference_rotation_boundary(
    boundary_reference::Matrix{Float64},
    boundary_target::Matrix{Float64},
    test_angles;
    metric=boundary_modified_hausdorff,
)
    shape_differences = Array{
        NamedTuple{(:angle, :shape_difference),Tuple{Float64,Float64}}
    }(
        undef, length(test_angles)
    )

    for (idx, angle) in enumerate(test_angles)
        # shape_difference_rotation undoes the candidate rotation with imrotate(-angle).
        # Boundary coordinates are image indices, whose row axis points down, so
        # rotate_boundary's positive (counterclockwise in an x-right, y-up frame) is the
        # opposite sense of imrotate's positive. Undoing imrotate(angle) is therefore
        # rotate_boundary(+angle). Verified against register on the same rotated mask.
        target_rotated = rotate_boundary(boundary_target, angle)
        shape_difference = metric(boundary_reference, target_rotated)
        shape_differences[idx] = (; angle, shape_difference)
    end
    return shape_differences
end

"""
    register_boundary(boundary_reference, boundary_target;
                      test_angles=register_default_angles_rad,
                      metric=boundary_modified_hausdorff)

Boundary-curve analogue of [`register`](@ref). Finds the angle in `test_angles` that
minimizes the shape difference between `boundary_reference` and `boundary_target`.

Shares `register`'s calling convention, so it can be passed straight to
`get_rotation_measurements(...; registration_function=register_boundary)` provided the
image column holds boundary matrices rather than masks. Ties resolve to whichever
candidate appears first in `test_angles`, matching `register`.
"""
function register_boundary(
    boundary_reference::Matrix{Float64},
    boundary_target::Matrix{Float64};
    test_angles=register_default_angles_rad,
    metric=boundary_modified_hausdorff,
)
    shape_differences = shape_difference_rotation_boundary(
        boundary_reference, boundary_target, test_angles; metric
    )
    best_match = argmin((x) -> x.shape_difference, shape_differences)
    return best_match.angle
end

"""
    boundary_shape_difference(boundary1, orientation1, boundary2, orientation2;
                              metric=boundary_modified_hausdorff)
    boundary_shape_difference(floe1::DataFrameRow, floe2::DataFrameRow; kwargs...)

Boundary-curve analogue of [`shape_difference`](@ref): aligns both curves on the same axis
using their estimated orientations, then compares them once with `metric`.

Like the mask-based version this performs a single comparison at the orientation-implied
pose rather than searching over angles — see [`register_boundary`](@ref) for the search.
That distinction matters for callers such as the candidate filters, which evaluate this
per candidate pair.

No centroid alignment step is needed: the metrics centre both curves internally, and
`x' - mean(x')` equals `R(x - mean(x))` regardless of the rotation centre, so translation
alignment is structural rather than something this function has to arrange.
"""
function boundary_shape_difference(
    boundary1::Matrix{Float64},
    orientation1::Real,
    boundary2::Matrix{Float64},
    orientation2::Real;
    metric=boundary_modified_hausdorff,
)
    # shape_difference applies imrotate(mask, orientation); in image coordinates that is
    # rotate_boundary(-orientation), see shape_difference_rotation_boundary.
    return metric(
        rotate_boundary(boundary1, -orientation1), rotate_boundary(boundary2, -orientation2)
    )
end

function boundary_shape_difference(
    boundary1::AbstractMatrix{<:Real},
    orientation1::Real,
    boundary2::AbstractMatrix{<:Real},
    orientation2::Real;
    kwargs...,
)
    return boundary_shape_difference(
        convert(Matrix{Float64}, boundary1),
        orientation1,
        convert(Matrix{Float64}, boundary2),
        orientation2;
        kwargs...,
    )
end

function boundary_shape_difference(floe1::DataFrameRow, floe2::DataFrameRow; kwargs...)
    return boundary_shape_difference(
        floe1.boundary, floe1.orientation, floe2.boundary, floe2.orientation; kwargs...
    )
end

"""
    BoundaryRegistration(; metric=boundary_modified_hausdorff, boundary_column=:boundary)

A configured, callable wrapper around [`register_boundary`](@ref) that fixes the shape
metric up front.

`get_rotation_measurements` invokes its `registration_function` as
`registration_function(image1, image2)`, with no way to thread a metric through, so a
metric that is not the default has to be baked into the callable instead. Because this
subtypes `Function` and its matrix method takes a `test_angles` keyword, an instance is a
drop-in `registration_function` whenever the image column holds boundary curves:

```julia
reg = BoundaryRegistration(; metric=boundary_hausdorff)
get_rotation_measurements(df; id_column=:id, image_column=:boundary,
                          time_column=:time, registration_function=reg)
```

Subtyping `Function` is required, not cosmetic: `get_rotation_measurements` annotates
`registration_function::Function`, and keyword arguments are converted, so a plain callable
struct would be rejected at the call site.

## Arguments
- `metric`: `metric(reference, rotated_target) -> Real`. See
  [`boundary_modified_hausdorff`](@ref) (default) and [`boundary_hausdorff`](@ref).
- `boundary_column`: column the `DataFrameRow` method reads.
"""
@kwdef struct BoundaryRegistration <: Function
    metric = boundary_modified_hausdorff
    boundary_column = :boundary
end

function (reg::BoundaryRegistration)(
    boundary_reference::Matrix{Float64},
    boundary_target::Matrix{Float64};
    test_angles=register_default_angles_rad,
)
    return register_boundary(
        boundary_reference, boundary_target; test_angles, metric=reg.metric
    )
end

# `:boundary` is an untyped column, so an Int matrix or a view can reach the functor. Widen
# rather than let it MethodError inside get_rotation_measurements' Threads.@threads, where
# it would surface as a TaskFailedException wrapping the real cause.
function (reg::BoundaryRegistration)(
    boundary_reference::AbstractMatrix{<:Real},
    boundary_target::AbstractMatrix{<:Real};
    kwargs...,
)
    return reg(
        convert(Matrix{Float64}, boundary_reference),
        convert(Matrix{Float64}, boundary_target);
        kwargs...,
    )
end

# Note: get_rotation_measurements never reaches this method -- it extracts
# row[image_column] first. This exists for direct use and for parity with
# `shape_difference(floe1::DataFrameRow, floe2::DataFrameRow)`.
function (reg::BoundaryRegistration)(floe1::DataFrameRow, floe2::DataFrameRow; kwargs...)
    return reg(floe1[reg.boundary_column], floe2[reg.boundary_column]; kwargs...)
end
