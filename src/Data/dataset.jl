export Case, Dataset, loader, metadata

import DataFrames: DataFrame, DataFrameRow, nrow, eachrow, subset
import FileIO: load, save

@kwdef struct Case
    loader::AbstractLoader
    metadata::DataFrameRow
end

function metadata(case::Case)::DataFrameRow
    return case.metadata
end

function loader(case::Case)::AbstractLoader
    return case.loader
end

struct Dataset
    loader::AbstractLoader
    metadata::DataFrame
end

function metadata(ds::Dataset)::DataFrame
    return ds.metadata
end

function loader(ds::Dataset)::AbstractLoader
    return ds.loader
end

Base.length(ds::Dataset)::Int = nrow(metadata(ds))
Base.iterate(ds::Dataset) = iterate(Case.(Ref(ds.loader), eachrow(metadata(ds))))
function Base.iterate(ds::Dataset, state)
    return iterate(Case.(Ref(loader(ds)), eachrow(metadata(ds))), state)
end
Base.filter(f::Function, ds::Dataset)::Dataset =
    Dataset(ds.loader, filter(row -> f(row), ds.metadata))
Base.getindex(ds::Dataset, i::Int)::Case = Case(ds.loader, ds.metadata[i, :])
subset(ds::Dataset, args...; kwargs...)::Dataset =
    Dataset(ds.loader, subset(ds.metadata, args...; kwargs...))