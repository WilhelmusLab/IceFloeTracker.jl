import Dates: Millisecond, Second, Minute, Hour, Day
abstract type AbstractTimeDistanceThresholdFunction end

# Minimum in Dates in milliseconds, so I just need to add
# a function that converts to base units.

function seconds(time_difference)
    ms = convert(Millisecond, time_difference)
    return ms.value / 1000
end

"""
LopezAcostaTimeDistanceFunction(Δx, Δt; dt, dx)

Stepwise time delta function based on Lopez-Acosta et al. 2019. The time thresholds and the input Δt
must be time objects so that conversion to seconds is possible. Displacement distances are assumed to
be in meters. The final dx value is the maximum displacement.
"""

@kwdef struct LopezAcostaTimeDistanceFunction <: AbstractTimeDistanceThresholdFunction
    dt=(Minute(20), Minute(90), Hour(24))
    dx=(3.75e3, 7.5e3, 30e3, 60e3)
    # TODO: add tests of input types and dimensions
end

function (f::LopezAcostaTimeDistanceFunction)(Δx, Δt)
    for (idx, time_threshold) in enumerate(f.dt)
        seconds(Δt) <= seconds(time_threshold) && Δx <= f.dx[idx] && return true
    end
    Δx <= f.dx[end] && return true
    return false
end

"""
Tests the travel distance and time in log-log space against an empirically fitted quadratic function. The
function is constrained by minimum and maximum times. Times less than the minimum are subject to the maximum 1-hour travel
distance, while times larger than the maximum fail automatically. See Watkins et al. 2025 for details.
"""

@kwdef struct LogLogQuadraticTimeDistanceFunction <: AbstractTimeDistanceThresholdFunction
    llq_params=[0.403, 0.988, -0.05]
    min_time=Hour(1)
    max_time=Day(7)
    # TODO: add tests of input types and dimensions
    # TODO: potentially move the parameterized equation into the struct
end

function (f::LogLogQuadraticTimeDistanceFunction)(Δx, Δt)
    a, b, c = f.llq_params
    seconds(Δt) > seconds(f.max_time) && return false
    
    Δx_km = Δx / 1000
    seconds(Δt) <= seconds(f.min_time) && return Δx_km <= 10^a # Check whether this is being activated early. Issue with time format?
    
    Δt_hours = seconds(Δt) / 3600
    Δx_km <= 10^(a + b * log10(Δt_hours) + c * (log10(Δt_hours))^2) && return true
    return false
end

"""
distance_threshold(Δx, Δt, threshold_function)

Time-distance threshold functions are used to identify search regions for floe matching. We include two distance threshold functions:
LopezAcostaTimeDistanceFunction, based on the stepwise method in Lopez-Acosta et al. 2019, and LogLogQuadraticTimeDistanceFunction, 
which is defined in Watkins et al. 2025 and is based on fitting a quadratic function to log-transformed displacements calculated from
drifting buoy data.

Example usage:
```
distance_threshold(100, Hour(12), LopezAcostaTimeDistanceFunction())
```

"""
# TODO: require dt to be milliseconds (or at least a timedelta), so we can do e.g. = Dates.seconds(passtimes[2] - passtimes[1])
function distance_threshold(Δx, Δt, threshold_function::AbstractTimeDistanceThresholdFunction)
    return threshold_function(Δx, Δt)
end
