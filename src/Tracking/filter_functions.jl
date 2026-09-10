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
julia> dt_test = DistanceThresholdFilter(time_colum=:Δt, dist_column=:Δx, threshold_function=LinearTimeDistanceFunction())
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
@kwdef struct DistanceThresholdFilter <: AbstractFloeFilterFunction
    time_column = :Δt
    dist_column = :Δx
    scaled_dist_column = :scaled_distance
    threshold_function = LinearTimeDistanceFunction()
    threshold_column = :time_distance_test
    scaling_function = dt -> maximum_linear_distance(dt; umax=1.5, eps=250)
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
    return sqrt.(
        (floe.row_centroid .- candidates.row_centroid) .^ 2 .+
        (floe.col_centroid .- candidates.col_centroid) .^ 2,
    ) * r
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
@kwdef struct RelativeErrorThresholdFilter <: AbstractFloeFilterFunction
    variable
    area_variable = :area
    threshold_column = :relative_error_test
    threshold_function = PiecewiseLinearThresholdFunction()
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
@kwdef struct ShapeDifferenceThresholdFilter <: AbstractFloeFilterFunction
    area_variable = :area
    scale_by = :area
    threshold_column = :shape_difference_test
    threshold_function = PiecewiseLinearThresholdFunction(100, 800, 0.5, 0.3)
end

function (f::ShapeDifferenceThresholdFilter)(
    floe::DataFrameRow, candidates::DataFrame, _::Val{:raw}
)
    function sd(mask, orientation)
        return round(
            shape_difference(floe.mask, floe.orientation, mask, orientation); digits=3
        )
    end

    transform!(candidates, [:mask, :orientation] => ByRow(sd) => :shape_difference)

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

"""
    PsiSCorrelationThresholdFunction(area_variable, threshold_column, threshold_function)
    PsiSCorrelationThresholdFunction(floe, candidates, Val(:raw))

Compute the ψ-s correlation between a floe and a dataframe of candidate floes. Adds the 
ψ-s correlation ``\\rho``,  ψ-s correlation score (1 - ``\\rho``), and the result of the threshold function
to the columns of `candidates`.
"""
@kwdef struct PsiSCorrelationThresholdFilter <: AbstractFloeFilterFunction
    area_variable = :area
    threshold_column = :psi_s_correlation_test
    threshold_function = PiecewiseLinearThresholdFunction(100, 800, 0.14, 0.1)
end

#TODO: Add option to include the confidence intervals with the normalized cross correlation tests.
function (f::PsiSCorrelationThresholdFilter)(floe, candidates, _::Val{:raw})
    rfloe(p2) = round(normalized_cross_correlation(floe.psi, p2); digits=3)
    transform!(candidates, [:psi] => ByRow(rfloe) => :psi_s_correlation)
    candidates[!, :psi_s_correlation_score] = 1 .- candidates[!, :psi_s_correlation]

    # Future work: add computation of the confidence intervals for psi-s corr here.
    return transform!(
        candidates,
        [f.area_variable, :psi_s_correlation_score] =>
            ByRow(f.threshold_function) => f.threshold_column,
    )
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

!!! warning "Thresholds are placeholders"
    The default 0.47 -> 0.31 values are inherited from the mask-based filter, where they
    were tuned for `count_symdiff / area` -- a ratio of pixel areas.
    `boundary_normalized_distance` is `MSE / perimeter^2`, a different dimensionless
    quantity on a much smaller scale, so these values admit nearly every candidate. They
    need re-tuning on real tracker data before this filter is used in anger, which is why
    it is not part of `default_filter`.

## Arguments
- `metric`: `metric(reference, rotated_target) -> Real`, default
  `boundary_normalized_distance`.
- `boundary_column`: boundary column to read, default `:boundary`.
- `area_variable`: column passed to `threshold_function` as the scale, default `:area`.

To let the score contribute to matching, opt into it alongside the filter:

```julia
cols = [MinimumWeightMatchingFunction().columns..., :scaled_boundary_shape_difference]
MinimumWeightMatchingFunction(; columns=cols, weights=ones(length(cols)))
```

`:scaled_boundary_shape_difference` is deliberately absent from that function's default
`columns`: a missing column there is not an error but a `@debug` log and an **empty**
result, so adding it while this filter stays opt-in would make the default tracker
silently match nothing.
"""
@kwdef struct BoundaryShapeDifferenceThresholdFilter <: AbstractFloeFilterFunction
    area_variable = :area
    boundary_column = :boundary
    metric = boundary_normalized_distance
    threshold_column = :boundary_shape_difference_test
    threshold_function = PiecewiseLinearThresholdFunction(100, 700, 0.47, 0.31)
end

function (f::BoundaryShapeDifferenceThresholdFilter)(
    floe::DataFrameRow, candidates::DataFrame, _::Val{:raw}
)
    reference = floe[f.boundary_column]
    function bsd(boundary, orientation)
        return boundary_shape_difference(
            reference, floe.orientation, boundary, orientation; metric=f.metric
        )
    end

    transform!(
        candidates,
        [f.boundary_column, :orientation] => ByRow(bsd) => :boundary_shape_difference,
    )

    # No division by :area here, deliberately. The mask filter scales because
    # count_symdiff returns a raw pixel count; boundary_normalized_distance is already
    # dimensionless, so scaling again would double-normalize. The column name is kept for
    # symmetry with the mask filter and with MinimumWeightMatchingFunction.
    candidates[!, :scaled_boundary_shape_difference] =
        candidates[!, :boundary_shape_difference]

    return transform!(
        candidates,
        [f.area_variable, :scaled_boundary_shape_difference] =>
            ByRow(f.threshold_function) => f.threshold_column,
    )
end

# Deliberately no `const boundary_shape_difference_filter` preset: filter_functions.jl is
# included before register.jl (Tracking.jl:49 vs :53), so constructing one at module scope
# would hit the `metric = boundary_normalized_distance` default before that function
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
