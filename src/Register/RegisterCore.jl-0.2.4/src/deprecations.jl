import Base: one
@deprecate one(::Type{NumDenom{T}}) where T     oneunit(NumDenom{T})
@deprecate one(p::NumDenom)                     oneunit(p)

@deprecate ratio(mm::AbstractArray, args...)  ratio.(mm, args...)
