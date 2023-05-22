module CenterIndexedArrays

using Interpolations, OffsetArrays
using OffsetArrays: IdentityUnitRange

export CenterIndexedArray

include("symrange.jl")

"""
A `CenterIndexedArray` is one for which the array center has indexes
`0,0,...`. Along each coordinate, allowed indexes range from `-n:n`.

CenterIndexedArray(A) "converts" `A` into a CenterIndexedArray. All
the sizes of `A` must be odd.
"""
struct CenterIndexedArray{T,N,A<:AbstractArray} <: AbstractArray{T,N}
    data::A
    halfsize::NTuple{N,Int}

    function CenterIndexedArray{T,N,A}(data::A) where {T,N,A<:AbstractArray}
        new{T,N,A}(data, _halfsize(data))
    end
end

CenterIndexedArray(A::AbstractArray{T,N}) where {T,N} = CenterIndexedArray{T,N,typeof(A)}(A)
CenterIndexedArray{T,N}(::UndefInitializer, sz::Vararg{<:Integer,N}) where {T,N} =
    CenterIndexedArray(Array{T,N}(undef, sz...))
CenterIndexedArray{T,N}(::UndefInitializer, sz::NTuple{N,Integer}) where {T,N} =
    CenterIndexedArray(Array{T,N}(undef, sz))
CenterIndexedArray{T}(::UndefInitializer, sz::Vararg{<:Integer,N}) where {T,N} =
    CenterIndexedArray{T,N}(undef, sz...)
CenterIndexedArray{T}(::UndefInitializer, sz::NTuple{N,Integer}) where {T,N} =
    CenterIndexedArray{T,N}(undef, sz)

# This is the AbstractArray default, but do this just to be sure
Base.IndexStyle(::Type{A}) where {A<:CenterIndexedArray} = IndexCartesian()

Base.size(A::CenterIndexedArray) = size(A.data)
Base.axes(A::CenterIndexedArray) = map(SymRange, A.halfsize)

const SymAx = Union{SymRange, Base.Slice{SymRange}}
Base.axes(r::Base.Slice{SymRange}) = (r.indices,)

function Base.similar(A::CenterIndexedArray, ::Type{T}, inds::Tuple{SymAx,Vararg{SymAx}}) where T
    data = Array{T}(undef, map(length, inds))
    CenterIndexedArray(data)
end
function Base.similar(::Type{T}, inds::Tuple{SymAx, Vararg{SymAx}}) where T
    data = Array{T}(undef, map(length, inds))
    CenterIndexedArray(data)
end

# This is incomplete: ideally we wouldn't need SymAx in the first slot
# as long as there was at least one SymAx.
function Base.similar(A::CenterIndexedArray, ::Type{T}, inds::Tuple{SymAx,Vararg{Union{Int,<:IdentityUnitRange,SymAx}}}) where T
    torange(n) = isa(n, Int) ? Base.OneTo(n) : n
    return OffsetArray{T}(undef, map(torange, inds))
end


function _halfsize(A::AbstractArray)
    all(isodd, size(A)) || error("Must have all-odd sizes")
    map(n->n>>UInt(1), size(A))
end

@inline function Base.getindex(A::CenterIndexedArray{T,N}, i::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, i...)
    @inbounds val = A.data[map(offset, A.halfsize, i)...]
    val
end

Base.@propagate_inbounds Base.getindex(A::CenterIndexedArray{T,N,I}, i::Vararg{Int,N}) where {T,N,I<:AbstractInterpolation} =
    _getindex(A, i...)
Base.@propagate_inbounds Base.getindex(A::CenterIndexedArray{T,N,I}, i::Vararg{Number,N}) where {T,N,I<:AbstractInterpolation} =
    _getindex(A, i...)

@inline function _getindex(A::CenterIndexedArray{T,N,I}, i::Vararg{Number,N}) where {T,N,I<:AbstractInterpolation}
    @boundscheck checkbounds(A, i...)
    @inbounds val = A.data(map(offset, A.halfsize, i)...)
    val
end
Base.throw_boundserror(A::CenterIndexedArray, I) = (Base.@_noinline_meta; throw(BoundsError(A, I)))

offset(off, i) = off+i+1

@inline function Base.setindex!(A::CenterIndexedArray{T,N}, v, i::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, i...)
    @inbounds A.data[map(offset, A.halfsize, i)...] = v
    v
end


Base.BroadcastStyle(::Type{<:CenterIndexedArray}) = Broadcast.ArrayStyle{CenterIndexedArray}()
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{CenterIndexedArray}}, ::Type{ElType}) where ElType
    similar(ElType, axes(bc))
end

Base.parent(A::CenterIndexedArray) = A.data

function Base.showarg(io::IO, A::CenterIndexedArray, toplevel)
    print(io, "CenterIndexedArray(")
    Base.showarg(io, parent(A), false)
    print(io, ')')
    toplevel && print(io, " with eltype ", eltype(A))
end

include("deprecated.jl")

end  # module
