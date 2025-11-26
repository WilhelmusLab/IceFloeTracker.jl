export Case, Dataset, loader, info

import DataFrames: DataFrames, DataFrame, DataFrameRow, nrow, eachrow, subset
import FileIO: load, save

@kwdef struct Case
    loader::AbstractLoader
    info::DataFrameRow
end

function info(case::Case)::DataFrameRow
    return case.info
end

function loader(case::Case)::AbstractLoader
    return case.loader
end

struct Dataset
    loader::AbstractLoader
    info::DataFrame
end

function info(ds::Dataset)::DataFrame
    return ds.info
end

function loader(ds::Dataset)::AbstractLoader
    return ds.loader
end

Base.length(ds::Dataset)::Int = nrow(info(ds))
Base.iterate(ds::Dataset) = iterate(Case.(Ref(ds.loader), eachrow(info(ds))))
function Base.iterate(ds::Dataset, state)
    return iterate(Case.(Ref(loader(ds)), eachrow(info(ds))), state)
end
Base.filter(f::Function, ds::Dataset)::Dataset =
    Dataset(ds.loader, filter(row -> f(row), ds.info))
Base.getindex(ds::Dataset, i::Int)::Case = Case(ds.loader, ds.info[i, :])
subset(ds::Dataset, args...; kwargs...)::Dataset =
    Dataset(ds.loader, subset(ds.info, args...; kwargs...))
