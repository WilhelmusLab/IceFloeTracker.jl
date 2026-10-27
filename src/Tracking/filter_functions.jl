"""
    AbstractFloeFilterFunction

The root type for the candidate filter functions.
"""
abstract type AbstractFloeFilterFunction <: Function end

function (f::AbstractFloeFilterFunction)(floe, candidates)
    f(floe, candidates, Val(:raw))
    subset!(candidates, [f.threshold_column] => r -> r .> 0)
    return select!(candidates, Not(f.threshold_column))
end

# TODO: Update the DistanceThresholdFilter to allow geospatial columns (i.e., call latlon first)
"""
    DistanceThresholdFilter(time_column, dist_column, threshold_function, threshold_column)
    DistanceThresholdFilter(floe, candidates)

The distance threshold filter creates columns for time and distance and applies a threshold
function to these columns to determine if the net travel is physically possible. The struct
is initialized with names for the time and distance columns, the threshold function (a TimeDistanceFunction)
and the name of the column in which to store the results.


```julia-repl
julia> dt_test = DistanceThresholdFilter(time_column=:Δt, dist_column=:Δx, threshold_function=LinearTimeDistanceFunction())
```
Now, let's assume that `floe` and `candidates` are already defined. Then

```julia-repl
julia> dt_test(floe, candidates)
```

will modify `candidates` in place to include only rows in which the `LinearTimeDistanceFunction()` evaluates as true.
Passing `Val{:raw}` as the third argument will forgo the subsetting step so that the output of the test can be examined.

## Arguments
- `time_column`: Name of the column to store pairwise floe time differences
- `dist_column`: Name of the column to store distances between floes
- `scaled_dist_column` = Name of the column to store scaled distances
- `threshold_function` = LinearTimeDistanceFunction()
- `threshold_column` = :time_distance_test
- `scaling_function`: Function producing a number (e.g., max expected distance) for scaling the distances. This quantify is used in the matching function as an error metric.
"""
@kwdef struct DistanceThresholdFilter{F,S} <: AbstractFloeFilterFunction
    time_column::Symbol = :Δt
    dist_column::Symbol = :Δx
    scaled_dist_column::Symbol = :scaled_distance
    threshold_function::F = LinearTimeDistanceFunction()
    threshold_column::Symbol = :time_distance_test
    scaling_function::S = dt -> maximum_linear_distance(dt; umax=1.5, eps=250)
end

function (f::DistanceThresholdFilter)(
    floe::DataFrameRow, candidates::DataFrame, _::Val{:raw}
) # can we get the same behavior with a less opaque function call?
    candidates[!, f.time_column] = candidates[!, :passtime] .- floe.passtime
    candidates[!, f.dist_column] = euclidean_distance(floe, candidates)
    candidates[!, f.scaled_dist_column] =
        candidates[!, f.dist_column] ./ f.scaling_function.(candidates[!, f.time_column])
    return transform!(
        candidates,
        [f.dist_column, f.time_column] => ByRow(f.threshold_function) => f.threshold_column,
    )
end

"""
    euclidean_distance(floe, candidates; r=250)

Compute the distance in meters between a floe and candidate floes by computing the
straight-line distance between centroids in pixel coordinates and converting that result
using a pixel resolution `r` with units meters/pixel. The floe and candidates must
have rows `row_centroid` and `col_centroid`.
"""
function euclidean_distance(floe, candidates; r=250)
    # Function barrier: the centroid columns are untyped, so read them once and compute
    # in a method that specializes on their types.
    return _euclidean_distance(
        floe.row_centroid,
        floe.col_centroid,
        candidates.row_centroid,
        candidates.col_centroid,
        r,
    )
end

function _euclidean_distance(row1, col1, rows2, cols2, r)
    return sqrt.((row1 .- rows2) .^ 2 .+ (col1 .- cols2) .^ 2) * r
end

# TODO: Add geodetic distance function

"""
    RelativeErrorThresholdFilter(variable, area_variable, threshold_column, threshold_function)
    RelativeErrorThresholdFilter(floe, candidates)
    RelativeErrorThresholdFilter(floe, candidates, Var(:raw))

Compute and test (absolute) relative error for `variable`. The relative error
between scalar variables X and Y is defined as
```math
\\eps = \\abs(X - Y)/\\text{mean}(X, Y)
```
This function takes a string or Symbol `variable` (which must be a named column in
the `candidates` DataFrame) and computes the relative error. Calling the function with
the variable name, `area_variable`, `threshold_column` name, and a `threshold_function`
initializes the function and saves the parameter values. Once initialized, the function
takes a `DataFrameRow` and a `DataFrame` of candidate floes as arguments, and subsets
the candidates to only those which evaluate as `true` using the `threshold_function`.
Including the dummy variable `Var(:raw)` returns the candidates dataframe with the test
results without subsetting it.
"""
@kwdef struct RelativeErrorThresholdFilter{V<:Union{Symbol,AbstractString},F} <:
              AbstractFloeFilterFunction
    variable::V
    area_variable::Symbol = :area
    threshold_column::Symbol = :relative_error_test
    threshold_function::F = PiecewiseLinearThresholdFunction()
end

function (f::RelativeErrorThresholdFilter)(
    floe::DataFrameRow, candidates::DataFrame, _::Val{:raw}
)
    new_variable = Symbol(:relative_error_, f.variable)
    X = floe[f.variable]
    Y = candidates[!, f.variable]
    candidates[!, new_variable] = abs.(X .- Y) ./ (0.5 .* (X .+ Y))
    return transform!(
        candidates,
        [f.area_variable, new_variable] =>
            ByRow(f.threshold_function) => f.threshold_column,
    )
end

"""
    ShapeDifferenceThresholdFilter(area_variable, scale_by, threshold_column, threshold_function)
    ShapeDifferenceThresholdFilter(floe, candidates)
    ShapeDifferenceThresholdFilter(floe, candidates, Val(:raw))


Compute and test the scaled shape difference between input `floe` and each floe in the dataframe `candidates`.
The shape difference between objects ``A`` and ``B`` is defined as
```math
SD = (A \\cup B) \\setminus (A \\cap B)
```
Here, the shapes are both rotated by their orientation and aligned at their respective centroids before computing ``\\SD``.
The result is divided by `scale_by` (e.g., area or perimeter), then the scaled value is assessed with
the `threshold_function` which is assumed to depend on area.

"""
@kwdef struct ShapeDifferenceThresholdFilter{F} <: AbstractFloeFilterFunction
    area_variable::Symbol = :area
    scale_by::Symbol = :area
    threshold_column::Symbol = :shape_difference_test
    threshold_function::F = PiecewiseLinearThresholdFunction(100, 800, 0.5, 0.3)
end

function (f::ShapeDifferenceThresholdFilter)(
    floe::DataFrameRow, candidates::DataFrame, _::Val{:raw}
)
    _add_shape_difference!(candidates, floe.mask, floe.orientation)

    candidates[!, :scaled_shape_difference] =
        candidates[!, :shape_difference] ./ candidates[!, f.scale_by]
    candidates[!, :scaled_shape_difference] .= round.(
        candidates[!, :scaled_shape_difference]; digits=3
    )

    return transform!(
        candidates,
        [f.area_variable, :scaled_shape_difference] =>
            ByRow(f.threshold_function) => f.threshold_column,
    )
end

# Function barrier: `floe` comes from untyped DataFrame columns, so its mask and orientation
# are read once by the caller and this method specializes on their types, instead of the
# per-candidate closure re-reading them from the row.
function _add_shape_difference!(candidates::DataFrame, mask1, orientation1)
    function sd(mask, orientation)
        return round(shape_difference(mask1, orientation1, mask, orientation); digits=3)
    end
    return transform!(candidates, [:mask, :orientation] => ByRow(sd) => :shape_difference)
end

"""
    PsiSCorrelationThresholdFunction(area_variable, threshold_column, threshold_function)
    PsiSCorrelationThresholdFunction(floe, candidates, Val(:raw))

Compute the ψ-s correlation between a floe and a dataframe of candidate floes. Adds the
ψ-s correlation ``\\rho``,  ψ-s correlation score (1 - ``\\rho``), and the result of the threshold function
to the columns of `candidates`.
"""
@kwdef struct PsiSCorrelationThresholdFilter{F} <: AbstractFloeFilterFunction
    area_variable::Symbol = :area
    threshold_column::Symbol = :psi_s_correlation_test
    threshold_function::F = PiecewiseLinearThresholdFunction(100, 800, 0.14, 0.1)
end

#TODO: Add option to include the confidence intervals with the normalized cross correlation tests.
function (f::PsiSCorrelationThresholdFilter)(floe, candidates, _::Val{:raw})
    _add_psi_s_correlation!(candidates, floe.psi)
    candidates[!, :psi_s_correlation_score] = 1 .- candidates[!, :psi_s_correlation]

    # Future work: add computation of the confidence intervals for psi-s corr here.
    return transform!(
        candidates,
        [f.area_variable, :psi_s_correlation_score] =>
            ByRow(f.threshold_function) => f.threshold_column,
    )
end

# Function barrier for the floe's ψ-s curve; see `_add_shape_difference!`.
function _add_psi_s_correlation!(candidates::AbstractDataFrame, psi1)
    rfloe(p2) = round(normalized_cross_correlation(psi1, p2); digits=3)
    return transform!(candidates, [:psi] => ByRow(rfloe) => :psi_s_correlation)
end

"""
    ChainedFilterFunction(filters::Vector{AbstractFloeFilterFunction})

A [`ChainedFilterFunction`](@ref) is a composite function based on a set of [`AbstractFloeFilterFunctions`](@ref). Each is
applied in sequence. Thus a filter function based on the distance threshold filter and area relative error filter
could be made as

```julia-repl
julia> filter_function = ChainedFilterFunction(
    filters=[DistanceThresholdFilter(), RelativeErrorThresholdFilter(variable=:area)]
)
```

Each item in the list should be an AbstractFloeFilterFunction.

"""
@kwdef struct ChainedFilterFunction <: AbstractFloeFilterFunction
    filters::Vector{AbstractFloeFilterFunction}
end

function (f::ChainedFilterFunction)(floe, candidates)
    for filter_fun in f.filters
        filter_fun(floe, candidates)
    end
end

"""
    BoundaryShapeDifferenceThresholdFilter(; kwargs...)

Boundary-curve counterpart to [`ShapeDifferenceThresholdFilter`](@ref). Compares each
candidate's boundary against the floe's using [`boundary_shape_difference`](@ref), which
aligns both by their orientations and applies `metric` once -- it does not search over
angles, matching the mask-based filter's cost model.

Requires a boundary column (see `add_boundary!`) and `:orientation` on both the floe and
the candidates. Writes `:boundary_shape_difference`,
`:scaled_boundary_shape_difference` and `threshold_column`.

## Arguments
- `threshold_function`: required. Called as `threshold_function(area, score)` with the
  score in pixels; returns whether the candidate passes.
- `metric`: `metric(reference, rotated_target) -> Real`, default
  `boundary_modified_hausdorff`.
- `boundary_column`: boundary column to read, default `:boundary`.
- `area_variable`: column passed to `threshold_function` as the scale, default `:area`.

`:scaled_boundary_shape_difference` is not among [`MinimumWeightMatchingFunction`](@ref)'s
default `columns`.
"""
@kwdef struct BoundaryShapeDifferenceThresholdFilter{M,F} <: AbstractFloeFilterFunction
    area_variable::Symbol = :area
    boundary_column::Symbol = :boundary
    metric::M = boundary_modified_hausdorff
    threshold_column::Symbol = :boundary_shape_difference_test
    # No default. ShapeDifferenceThresholdFilter's bounds (0.47 -> 0.31) were calibrated
    # against count_symdiff / area, a dimensionless ratio of pixel areas; the score here is
    # a distance in PIXELS, so they are wrong in units, not only in scale. Bounds in pixels,
    # area-dependent via PiecewiseLinearThresholdFunction, need to be estimated on tracked
    # data. Until then the filter is excluded from default_filter.
    threshold_function::F
end

function (f::BoundaryShapeDifferenceThresholdFilter)(
    floe::DataFrameRow, candidates::DataFrame, _::Val{:raw}
)
    _add_boundary_shape_difference!(
        candidates, f.boundary_column, floe[f.boundary_column], floe.orientation, f.metric
    )

    # No division by :area here. The mask filter scales because count_symdiff returns a
    # raw pixel count; a Hausdorff-type score is already a length, and how "similar" should
    # depend on floe size belongs in threshold_function, which takes area as an input. The
    # column name is kept for symmetry with the mask filter and MinimumWeightMatchingFunction.
    # TODO: not used for matching. The score is in pixels, while
    # MinimumWeightMatchingFunction sums its columns unweighted and the others are
    # dimensionless; whether and how to include it in the ranking is untested.
    candidates[!, :scaled_boundary_shape_difference] =
        candidates[!, :boundary_shape_difference]

    return transform!(
        candidates,
        [f.area_variable, :scaled_boundary_shape_difference] =>
            ByRow(f.threshold_function) => f.threshold_column,
    )
end

# Function barrier for the floe's boundary and orientation; see `_add_shape_difference!`.
function _add_boundary_shape_difference!(
    candidates::DataFrame, boundary_column::Symbol, reference, orientation1, metric::M
) where {M} # `where` forces specialization on a function-valued `metric`
    function bsd(boundary, orientation)
        return boundary_shape_difference(
            reference, orientation1, boundary, orientation; metric
        )
    end
    return transform!(
        candidates,
        [boundary_column, :orientation] => ByRow(bsd) => :boundary_shape_difference,
    )
end

# Deliberately no `const boundary_shape_difference_filter` preset: filter_functions.jl is
# included before register.jl (Tracking.jl:49 vs :53), so constructing one at module scope
# would hit the `metric = boundary_modified_hausdorff` default before that function
# exists. @kwdef evaluates defaults at construction time, so the struct itself is fine.

const max_travel_distance_filter = DistanceThresholdFilter(;
    threshold_function=LogLogQuadraticTimeDistanceFunction()
)

const area_relative_error_filter = RelativeErrorThresholdFilter(;
    variable=:area,
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.43, maximum_value=0.17
    ),
)

const convex_area_relative_error_filter = RelativeErrorThresholdFilter(;
    variable=:convex_area,
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.44, maximum_value=0.25
    ),
)

const major_axis_relative_error_filter = RelativeErrorThresholdFilter(;
    variable=:major_axis_length,
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.27, maximum_value=0.13
    ),
)

const minor_axis_relative_error_filter = RelativeErrorThresholdFilter(;
    variable=:minor_axis_length,
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.28, maximum_value=0.1
    ),
)

const shape_difference_filter = ShapeDifferenceThresholdFilter(;
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.47, maximum_value=0.31
    ),
)

const psi_s_correlation_filter = PsiSCorrelationThresholdFilter(;
    threshold_function=PiecewiseLinearThresholdFunction(;
        minimum_area=100, maximum_area=700, minimum_value=0.86, maximum_value=0.96
    ),
)

const default_filter = [
    max_travel_distance_filter,
    area_relative_error_filter,
    convex_area_relative_error_filter,
    major_axis_relative_error_filter,
    minor_axis_relative_error_filter,
    shape_difference_filter,
    psi_s_correlation_filter,
]

"""FilterFunction()

The default filter function for the FloeTracker. The function is an instance of [`ChainedFilterFunction`](@ref),
applying 7 individual [`AbstractFloeFilterFunctions`](@ref) in sequence:
    1. `DistanceThresholdFilter`
    2. `RelativeErrorThresholdFilters` for area, convex area, major axis length, and minor axis length
    3. `ShapeDifferenceThresholdFilter`
    4. `PsiSCorrelationThresholdFilter`
Filters in step 2 use the [`PiecewiseLinearThresholdFunction`](@ref) for thresholds, while Filter 1 uses a [`LinearTimeDistanceFunction`](@ref).
""" # TODO: Add reference to cal-val paper when ready.
function FilterFunction()
    ChainedFilterFunction(; filters=default_filter)
end

"""
    LopezAcosta2019ChainedFilterFunction(floe, candidates)

The LopezAcosta2019ChainedFilterFunction is a special case of ChainedFilterFunction
with parameters and threshold functions set based on Lopez-Acosta et al. 2019. The set
of threshold filters is the same as in the default FilterFunction, but using
stepwise threshold functions instead of piecewise.
"""
LopezAcosta2019ChainedFilterFunction = ChainedFilterFunction(;
    filters=[
        DistanceThresholdFilter(; threshold_function=LopezAcostaTimeDistanceFunction()),
        RelativeErrorThresholdFilter(; variable=:area), # use the step functions
        RelativeErrorThresholdFilter(; variable=:convex_area),
        RelativeErrorThresholdFilter(; variable=:major_axis_length),
        RelativeErrorThresholdFilter(; variable=:minor_axis_length),
        ShapeDifferenceThresholdFilter(), # Replace with step function
        PsiSCorrelationThresholdFilter(), # Replace with step function
    ],
)
