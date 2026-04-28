using DataFrames
using Dates
using Tables

candidates = DataFrame(;
    head_uuid=["f1", "f1", "f2", "f2"],
    uuid=["g1", "g2", "g1", "g2"],
    passtime=[
        DateTime(2020, 1, 1, 12, 0),
        DateTime(2020, 1, 1, 12, 5),
        DateTime(2020, 1, 1, 12, 10),
        DateTime(2020, 1, 1, 12, 15),
    ],
);
floe = Tables.Row((; :passtime => DateTime(2020, 1, 1, 12, 0)));

## Simplest case: no special parameters, just a function of time difference
function ΔPasstimeTransform(floe)
    return :passtime => ByRow((t) -> abs(t - floe.passtime)) => :Δpasstime
end

function ΔPasstimeHalfTransform(floe)
    return :passtime => ByRow((t) -> -(t - floe.passtime) / 2) => :Δpasstime∇
end

@show ΔPasstimeTransform(floe)
@show transform(candidates, ΔPasstimeTransform(floe), ΔPasstimeHalfTransform(floe))

## Next case: a parametrized transform function, where the parameter is set as a default
function ΔPasstimeTransformEps(floe; eps=Hour(6))
    return :passtime => ByRow((t) -> abs(t - floe.passtime) + eps) => :Δpasstime
end

@show ΔPasstimeTransformEps(floe)
@show transform(candidates, ΔPasstimeTransformEps(floe))

## Next case: a parametrized transform function, where the parameter is set via a functor struct

@kwdef struct ΔPasstimeTransformEpsFunctor
    eps::Period = Hour(6)
end

function (f::ΔPasstimeTransformEpsFunctor)(floe)
    return :passtime => ByRow((t) -> abs(t - floe.passtime) + f.eps) => :Δpasstime
end

@show ΔPasstimeTransformEpsFunctor()(floe)
@show transform(candidates, ΔPasstimeTransformEpsFunctor()(floe))

## Now add a filtering step based on the transform

ΔPasstimeFilter = [:Δpasstime] => ByRow((Δt) -> Δt < Hour(16))
@show dftf = transform(candidates, ΔPasstimeTransformEpsFunctor()(floe))
@show subset(dftf, ΔPasstimeFilter)
