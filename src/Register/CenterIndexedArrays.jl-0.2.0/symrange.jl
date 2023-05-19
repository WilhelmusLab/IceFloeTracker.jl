## SymRange, an AbstractUnitRange that's symmetric around 0
# These are used as axes for CenterIndexedArrays
struct SymRange <: AbstractUnitRange{Int}
    n::Int  # goes from -n:n
end

function SymRange(r::AbstractUnitRange)
    first(r) == -last(r) || error("cannot convert $r to a SymRange")
    SymRange(last(r))
end

Base.first(r::SymRange) = -r.n
Base.last(r::SymRange) = r.n
Base.axes(r::SymRange) = (r,)

@inline Base.unsafe_indices(r::SymRange) = (r,)

function iterate(r::SymRange)
    r.n == 0 && return nothing
    first(r), first(r)
end

function iterate(r::SymRange, s)
    s == last(r) && return nothing
    copy(s+1), s+1
end

@inline function Base.getindex(v::SymRange, i::Int)
    @boundscheck abs(i) <= v.n || Base.throw_boundserror(v, i)
    return i
end

Base.intersect(r::SymRange, s::SymRange) = SymRange(min(last(r), last(s)))

@inline function Base.getindex(r::SymRange, s::SymRange)
    @boundscheck checkbounds(r, s)
    s
end

@inline function Base.getindex(r::SymRange, s::AbstractUnitRange{<:Integer})
    @boundscheck checkbounds(r, s)
    return s
end

# TODO: should we be worried about the mismatch in axes?
# And should `convert(SymRange, r)` fail if axes(r) isn't the same as the result?
Base.promote_rule(::Type{SymRange}, ::Type{UR}) where {UR<:AbstractUnitRange} =
    UR
Base.promote_rule(::Type{UnitRange{T2}}, ::Type{SymRange}) where {T2} =
    UnitRange{promote_type(T2, Int)}

Base.show(io::IO, r::SymRange) = print(io, "SymRange(", repr(last(r)), ')')

if isdefined(Base, :reduced_index)
    Base.reduced_index(r::SymRange) = SymRange(0)
end
