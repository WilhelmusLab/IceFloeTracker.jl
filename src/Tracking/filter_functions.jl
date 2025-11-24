# TODO
# Set default functions for all filter functions based on the calibration paper results

"""
    AbstractFloeFilterFunction

The root type for the candidate filter functions.
"""
abstract type AbstractFloeFilterFunction <: Function end

function (f::AbstractFloeFilterFunction)(floe, candidates)
    f(floe, candidates, Val(:raw))
    subset!(candidates, [f.threshold_column] => r -> r .> 0)
    select!(candidates, Not(f.threshold_column))
end

# TODO: Update the DistanceThresholdFilter to allow geospatial columns (i.e., call latlon first)
"""
    DistanceThresholdFilter(time_column, dist_column, threshold_function, threshold_column)
    DistanceThresholdFilter(floe, candidates)

The distance threshold filter creates columns for time and distance and applies a threshold
function to these columns to determine if the net travel is physically possible. The struct
is initialized with names for the time and distance columns, the threshold function (a TimeDistanceFunction)
and the name of the column in which to store the results. 


```
julia> dt_test = DistanceThresholdFilter(time_colum=:Δt, dist_column=:Δx, threshold_function=LinearTimeDistanceFunction())
```
Now, let's assume that `floe` and `candidates` are already defined. Then

```
julia> dt_test(floe, candidates)
```

will modify `candidates` in place to include only rows in which the `LinearTimeDistanceFunction()` evaluates as true. 
Passing `Val{:raw}` as the third argument will forgo the subsetting step so that the output of the test can be examined.

""" 
@kwdef struct DistanceThresholdFilter <: AbstractFloeFilterFunction
        time_column = :Δt
        dist_column = :Δx
        scaled_dist_column = :scaled_distance
        threshold_function = LinearTimeDistanceFunction()
        threshold_column = :time_distance_test
        scaling_function = dt -> maximum_linear_distance(dt; umax=1.5, eps=250)
end

function (f::DistanceThresholdFilter)(floe::DataFrameRow, candidates::DataFrame, _::Val{:raw}) # can we get the same behavior with a less opaque function call?
    candidates[!, f.time_column] = candidates[!, :passtime] .- floe.passtime
    candidates[!, f.dist_column] = euclidean_distance(floe, candidates)
    candidates[!, f.scaled_dist_column] = candidates[!, f.dist_column] ./ f.scaling_function.(candidates[!, f.time_column])
    transform!(candidates, [f.dist_column, f.time_column] => 
        ByRow(f.threshold_function) => f.threshold_column)
end


"""
    euclidean_distance(floe, candidates; r=250)

Compute the distance in meters between a floe and candidate floes by computing the
straight-line distance between centroids in pixel coordinates and converting that result
using a pixel resolution `r` with units meters/pixel.
"""
function euclidean_distance(floe, candidates; r = 250)
    return sqrt.((floe.row_centroid .- candidates.row_centroid).^2 .+ 
    (floe.col_centroid .- candidates.col_centroid).^2) * r
end


"""
    RelativeErrorThresholdFilter(variable, area_variable, threshold_column, threshold_function)
    RelativeErrorThresholdFilter(floe, candidates)
    RelativeErrorThresholdFilter(floe, candidates, Var(:raw))

Compute and test (absolute) relative error for `variable`. The relative error
between scalar variables X and Y is defined as 
```
err = abs(X - Y)/mean(X, Y)
```
This function takes a scalar `<variable>` (which must be a named column in 
the `candidates` DataFrame) and computes the relative error. Calling the function with 
the variable name, `area_variable`, `threshold_column name``, and a `threshold_function`
initializes the function and saves the parameter values. Once initialized, the function 
takes a floe (DataFrameRow) and a DataFrame of candidate floes as arguments, and subsets
the candidates to only those which evaluate as true using the `threshold_function`.
Including the dummy variable `Var(:raw)` returns the candidates dataframe with the test 
results without subsetting it.
"""
@kwdef struct RelativeErrorThresholdFilter <: AbstractFloeFilterFunction
    variable
    area_variable = :area
    threshold_column = :relative_error_test
    threshold_function = PiecewiseLinearThresholdFunction()
end

function (f::RelativeErrorThresholdFilter)(floe::DataFrameRow, candidates::DataFrame, _::Val{:raw})
    new_variable = Symbol(:relative_error_, f.variable)
    X = floe[f.variable]
    Y = candidates[!, f.variable]
    candidates[!, new_variable] = abs.(X .- Y) ./ (0.5 .* (X .+ Y))
    transform!(candidates, [f.area_variable, new_variable] => 
        ByRow(f.threshold_function) => f.threshold_column)
end

"""
    ShapeDifferenceThresholdFilter(area_variable, scale_by, threshold_column, threshold_function)
    ShapeDifferenceThresholdFilter(floe, candidates)
    ShapeDifferenceThresholdFilter(floe, candidates, Val(:raw))
    

Compute and test the scaled shape difference between input `floe` and each floe in the dataframe `candidates`.
Assumes that the shape difference test operates on the shape difference scaled by a variable `scale_by`
and the shape difference test depends on the area. 

"""
@kwdef struct ShapeDifferenceThresholdFilter <: AbstractFloeFilterFunction
    area_variable = :area
    scale_by = :area
    threshold_column = :shape_difference_test
    threshold_function = PiecewiseLinearThresholdFunction(100, 800, 0.5, 0.3)
end

function (f::ShapeDifferenceThresholdFilter)(floe::DataFrameRow, candidates::DataFrame, _::Val{:raw})
    sd(mask, orientation) = round(shape_difference(floe.mask, floe.orientation, mask, orientation), digits=3)
    
    transform!(candidates,  [:mask, :orientation] => 
        ByRow(sd) => :shape_difference)

    candidates[!, :scaled_shape_difference] = candidates[!, :shape_difference] ./ candidates[!, f.scale_by]
    candidates[!, :scaled_shape_difference] .= round.(candidates[!, :scaled_shape_difference], digits=3)
    
    transform!(candidates, [f.area_variable, :scaled_shape_difference] =>
        ByRow(f.threshold_function) => f.threshold_column)
end

"""
    PsiSCorrelationThresholdFunction(area_variable, threshold_column, threshold_function)
    PsiSCorrelationThresholdFunction(floe, candidates, Val(:raw))

Compute the psi-s correlation between a floe and a dataframe of candidate floes. Adds the 
psi-s correlation,  psi-s correlation score (1 - correlation), and the result of the threshold function
to the columns of `candidates`.
""" 
@kwdef struct PsiSCorrelationThresholdFilter <: AbstractFloeFilterFunction
    area_variable = :area
    threshold_column = :psi_s_correlation_test
    threshold_function = PiecewiseLinearThresholdFunction(100, 800, 0.14, 0.1)
end

#TODO: Add option to include the confidence intervals with the normalized cross correlation tests.
function (f::PsiSCorrelationThresholdFilter)(floe, candidates, _::Val{:raw})
    

    rfloe(p2) = round(normalized_cross_correlation(floe.psi, p2), digits=3)
    transform!(candidates,  [:psi] => ByRow(rfloe) => :psi_s_correlation)
    candidates[!, :psi_s_correlation_score] = 1 .- candidates[!, :psi_s_correlation]

    # Future work: add computation of the confidence intervals for psi-s corr here.
    transform!(candidates, [f.area_variable, :psi_s_correlation_score] =>
        ByRow(f.threshold_function) => f.threshold_column
    )
end

"""
    ChainedFilterFunction(filters::Vector{AbstractFloeFilterFunction})

A ChainedFilterFunction is a composite function based on a set of AbstractFloeFilterFunctions. Each is
applied in sequence. Thus a filter function based on the DistanceThresholdFilter and area relative error 
could be made as

```julia
filter_function = ChainedFilterFunction(
    filters=[DistanceThresholdFilter(), RelativeErrorThresholdFilter(variable=:area)]
)
```

""" 
@kwdef struct ChainedFilterFunction <: AbstractFloeFilterFunction
    filters::Vector{AbstractFloeFilterFunction}
end

function (f::ChainedFilterFunction)(floe, candidates)
    for filter_fun in f.filters
        filter_fun(floe, candidates)
    end
end

"""FilterFunction()

The default filter function for the FloeTracker. The function is an instance of ChainedFilterFunction,
applying 7 individual AbstractFloeFilterFunctions in sequence:
    1. DistanceThresholdFilter
    2-5. RelativeErrorThresholdFilters for area, convex_area, major_axis_length, and minor_axis_length
    6. ShapeDifferenceThresholdFilter
    7. PsiSCorrelationThresholdFilter
Filters 2-7 use PiecewiseLinearThresholdFunctions for thresholds, while Filter 1 uses a LinearTimeDistanceFunction.
The default values and settings are derived in Watkins et al. 2026.
"""
FilterFunction = ChainedFilterFunction(
    filters = [
        DistanceThresholdFilter(),
        RelativeErrorThresholdFilter(variable=:area,
             threshold_function=PiecewiseLinearThresholdFunction(minimum_value=0.43, maximum_value=0.17)),
        RelativeErrorThresholdFilter(variable=:convex_area, 
            threshold_function=PiecewiseLinearThresholdFunction(minimum_value=0.44, maximum_value=0.25)),
        RelativeErrorThresholdFilter(variable=:major_axis_length, 
            threshold_function=PiecewiseLinearThresholdFunction(minimum_value=0.27, maximum_value=0.13)),
        RelativeErrorThresholdFilter(variable=:minor_axis_length, 
            threshold_function=PiecewiseLinearThresholdFunction(minimum_value=0.28, maximum_value=0.1)),
        ShapeDifferenceThresholdFilter(
            threshold_function=PiecewiseLinearThresholdFunction(minimum_value=0.47, maximum_value=0.31)),
        PsiSCorrelationThresholdFilter(
            threshold_function=PiecewiseLinearThresholdFunction(minimum_value=0.86, maximum_value=0.96))
        ])

"""
    LopezAcosta2019ChainedFilterFunction(floe, candidates)

The LopezAcosta2019ChainedFilterFunction is a special case of ChainedFilterFunction
with parameters and threshold functions set based on Lopez-Acosta et al. 2019. The set
of threshold filters is the same as in the default FilterFunction, but using
stepwise threshold functions instead of piecewise.
"""
LopezAcosta2019ChainedFilterFunction = ChainedFilterFunction(
    filters=[
        DistanceThresholdFilter(threshold_function=LopezAcostaTimeDistanceFunction()),
        RelativeErrorThresholdFilter(variable=:area), # use the step functions
        RelativeErrorThresholdFilter(variable=:convex_area),
        RelativeErrorThresholdFilter(variable=:major_axis_length),
        RelativeErrorThresholdFilter(variable=:minor_axis_length),
        ShapeDifferenceThresholdFilter(),
        PsiSCorrelationThresholdFilter()
        ]
)