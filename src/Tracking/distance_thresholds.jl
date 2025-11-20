import Dates: Millisecond, Second, Minute, Hour, Day
abstract type AbstractThresholdFunction <: Function end
abstract type AbstractTimeDistanceThresholdFunction <: AbstractThresholdFunction end

# Minimum in Dates in milliseconds, so I just need to add
# a function that converts to base units.

"""
    seconds(time_difference)

Convenience function to convert time difference in Millisecond to decimal seconds.
"""
function seconds(time_difference::Union{Period,CompoundPeriod})
    typeof(time_difference) == Dates.Millisecond && return time_difference.value / 1000
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
Tests the travel distance and time to a linear estimate of maximum travel distance
using the formula
```
max_dx = max_vel * dt + eps
```
Epsilon should be the uncertainty in position, such that if for example the positional
uncertainty is 250 m, then the maximum distance includes a 250 m buffer. The default maximum
velocity is 1.5 m/s.

"""

@kwdef struct LinearTimeDistanceFunction <: AbstractTimeDistanceThresholdFunction
    max_velocity = 0.75
    epsilon = 250
end

function (f::LinearTimeDistanceFunction)(Δx, Δt)
    max_Δx = maximum_linear_distance(Δt; umax=f.max_velocity, eps=f.epsilon)
    return max_Δx > Δx
end

"""
    maximum_linear_distance(Δt; umax=2, eps=250)

Compute the maximum travel distance based on travel time Δt (a Time Period) based on
the maximum velocity `umax` and an additive uncertainty of `eps` meters.
"""
function maximum_linear_distance(Δt; umax=0.75, eps=250)
    s = seconds(Δt)
    return s*umax + eps
end


"""
distance_threshold(Δx, Δt, threshold_function)

Time-distance threshold functions are used to identify search regions for floe matching. We include two distance threshold functions:
LopezAcostaTimeDistanceFunction, based on the stepwise method in Lopez-Acosta et al. 2019, and LogLogQuadraticTimeDistanceFunction, 
which is based on fitting a quadratic function to log-transformed displacements calculated from
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
