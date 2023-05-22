@deprecate CenterIndexedArray(::Type{T}, dims) where {T}      CenterIndexedArray{T}(undef, dims...)
@deprecate CenterIndexedArray(::Type{T}, dims...) where {T}   CenterIndexedArray{T}(undef, dims...)
