# TODO: Update to take square root of area (parameterized in terms of length scale)

"""
The piecewise linear threshold function is defined using two (area, value) pairs. For
areas below the minimum area, it is constant at minimum value; likewise for above the
maximum area. The threshold function is linear in between these two points. A return 
value `true` indicates that the value is below the threshold. 
"""
@kwdef struct PiecewiseLinearThresholdFunction <: AbstractThresholdFunction
    minimum_area = 100
    maximum_area = 700
    minimum_value = 0.4
    maximum_value = 0.2
end

function (f::PiecewiseLinearThresholdFunction)(area, value)
    area < f.minimum_area && return value < f.minimum_value
    area > f.maximum_area && return value < f.maximum_value
    slope = (f.maximum_value - f.minimum_value) / (f.maximum_area - f.minimum_area)
    return value < slope*(area - f.maximum_area) + f.maximum_value
end

"""
The stepwise linear threshold function is defined using a changepoint area and two levels. 
If the area is less than the changepoint area, the function returns true if the value is below
`low_value` and false otherwsie; if the area is greater than or equal to the changepoint area, 
then the value is tested againg `high_value`.
"""
@kwdef struct StepwiseLinearThresholdFunction <: AbstractThresholdFunction
    changepoint_area::Number
    low_value::Number
    high_value::Number
end

function (f::StepwiseLinearThresholdFunction)(area, value)
    area < f.changepoint_area && return value < f.low_value
    area >= f.changepoint_area && return value < f.high_value
end
