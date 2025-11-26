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

Base.getindex(ds::Dataset, i::Int) = Case(ds.loader, ds.info[i, :])

Base.iterate(ds::Dataset) = iterate(Case.(Ref(ds.loader), eachrow(info(ds))))
function Base.iterate(ds::Dataset, state)
    return iterate(Case.(Ref(loader(ds)), eachrow(info(ds))), state)
end

Base.filter(f::Function, ds::Dataset) = Dataset(ds.loader, filter(row -> f(row), ds.info))
Base.filter!(f::Function, ds::Dataset) = filter!(row -> f(row), ds.info)

function subset(ds::Dataset, args...; kwargs...)
    return Dataset(ds.loader, subset(ds.info, args...; kwargs...))
end
subset!(ds::Dataset, args...; kwargs...) = subset!(ds.info, args...; kwargs...)

function Base.sort(ds::Dataset, args...; kwargs...)
    return Dataset(ds.loader, sort(ds.info, args...; kwargs...))
end
function Base.sort(ds::Dataset, args...; kwargs...)
    return Dataset(ds.loader, sort(ds.info, args...; kwargs...))
end
function Base.sort(cols::Vector{Symbol}, ds::Dataset; kwargs...)
    return Dataset(ds.loader, sort(cols, ds.info; kwargs...))
end
Base.sort!(cols::Vector{Symbol}, ds::Dataset; kwargs...) = sort!(ds.info, cols; kwargs...)
Base.sort!(ds::Dataset, args...; kwargs...) = sort!(ds.info, args...; kwargs...)

Base.reverse(ds::Dataset) = Dataset(ds.loader, reverse(ds.info))
Base.reverse!(ds::Dataset) = reverse!(ds.info)

