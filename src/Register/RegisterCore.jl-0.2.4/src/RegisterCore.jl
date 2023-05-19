module RegisterCore

using CenterIndexedArrays
using ImageCore, ImageFiltering
using Requires

import Base: +, -, *, /

export
    # types
    MismatchArray,
    NumDenom,
    ColonFun,
    PreprocessSNF,
    # functions
    highpass,
    indmin_mismatch,
    maxshift,
    mismatcharrays,
    ratio,
    separate,
    paddedview,
    trimmedview

"""
The major functions/types exported by this module are:

- `NumDenom` and `MismatchArray`: packed pair representation of
  `(num,denom)` mismatch data
- `separate`: splits a `NumDenom` array into its component `num,denom` arrays
- `indmin_mismatch`: find the location of the minimum mismatch
- `highpass`: highpass filter an image before performing registration
- `PreprocessSNF`: shot-noise standardization and filtering

"""
RegisterCore

"""
    x = NumDenom(num, denom)

`NumDenom{T}` is an object containing a `(num,denom)` pair.
`x.num` is `num` and `x.denom` is `denom`.

Algebraically, `NumDenom` objects act like 2-vectors, and can be added and multiplied
by scalars:

    nd1 + nd2 = NumDenom(nd1.num + nd2.num, nd1.denom + nd2.denom)
    2*nd      = NumDenom(2*nd.num, 2*nd.denom)

Note that this is *not* what you'd get from normal arithmetic with ratios, where, e.g.,
`2*nd` would be expected to produce `NumDenom(2*nd.num, nd.denom)`.
The reason for calling these `*` and `+` is for use in `Interpolations.jl`, because it
allows interpolation to be performed on "both arrays" at once without
recomputing the interpolation coefficients. See the documentation for information
about how this is used for performing aperturered mismatch computations.

As a consequence, there is no `convert(Float64, nd::NumDenom)` method, because the algebra
above breaks any pretense that `NumDenom` numbers are somehow equivalent to ratios.
If you want to convert to a ratio, see [`ratio`](@ref).
"""
struct NumDenom{T<:Number}
    num::T
    denom::T
end
NumDenom(n::Gray, d::Gray) = NumDenom(gray(n), gray(d))
NumDenom(n::Gray, d) = NumDenom(gray(n), d)
NumDenom(n, d::Gray) = NumDenom(n, gray(d))
NumDenom(n, d) = NumDenom(promote(n, d)...)

(+)(p1::NumDenom, p2::NumDenom) = NumDenom(p1.num+p2.num, p1.denom+p2.denom)
(-)(p1::NumDenom, p2::NumDenom) = NumDenom(p1.num-p2.num, p1.denom-p2.denom)
(*)(n::Number, p::NumDenom) = NumDenom(n*p.num, n*p.denom)
(*)(p::NumDenom, n::Number) = n*p
(/)(p::NumDenom, n::Number) = NumDenom(p.num/n, p.denom/n)
Base.oneunit(::Type{NumDenom{T}}) where {T} = NumDenom{T}(oneunit(T),oneunit(T))
Base.oneunit(p::NumDenom) = oneunit(typeof(p))
Base.zero(::Type{NumDenom{T}}) where {T} = NumDenom(zero(T),zero(T))
Base.zero(p::NumDenom) = zero(typeof(p))
Base.promote_rule(::Type{NumDenom{T1}}, ::Type{T2}) where {T1,T2<:Number} = NumDenom{promote_type(T1,T2)}
Base.promote_rule(::Type{NumDenom{T1}}, ::Type{NumDenom{T2}}) where {T1,T2} = NumDenom{promote_type(T1,T2)}
Base.eltype(::Type{NumDenom{T}}) where {T} = T
Base.convert(::Type{NumDenom{T}}, p::NumDenom{T}) where {T} = p
Base.convert(::Type{NumDenom{T}}, p::NumDenom) where {T} = NumDenom{T}(p.num, p.denom)

Base.round(::Type{NumDenom{T}}, p::NumDenom) where T = NumDenom{T}(round(T, p.num), round(T, p.denom))

Base.convert(::Type{F}, p::NumDenom) where F<:AbstractFloat = error("`convert($F, ::NumDenom)` is deliberately not defined, see `?NumDenom`.")

function Base.show(io::IO, p::NumDenom)
    print(io, "NumDenom(")
    show(io, p.num)
    print(io, ",")
    show(io, p.denom)
    print(io, ")")
    return
end

const MismatchArray{ND<:NumDenom,N,A} = CenterIndexedArray{ND,N,A}

"""
    mxs = maxshift(D)

Return the `maxshift` value used to compute the mismatch array `D`.
"""
maxshift(A::MismatchArray) = A.halfsize

"""
`numdenom = MismatchArray(num, denom)` packs the array-pair
`(num,denom)` into a single `MismatchArray`.  This is useful
preparation for interpolation.
"""
function (::Type{M})(num::AbstractArray, denom::AbstractArray) where M<:MismatchArray
    size(num) == size(denom) || throw(DimensionMismatch("num and denom must have the same size"))
    T = promote_type(eltype(num), eltype(denom))
    numdenom = CenterIndexedArray{NumDenom{T}}(undef, size(num))
    _packnd!(numdenom, num, denom)
end

function _packnd!(numdenom::AbstractArray, num::AbstractArray, denom::AbstractArray)
    Rnd, Rnum, Rdenom = eachindex(numdenom), eachindex(num), eachindex(denom)
    if Rnum == Rdenom
        for (Idest, Isrc) in zip(Rnd, Rnum)
            @inbounds numdenom[Idest] = NumDenom(num[Isrc], denom[Isrc])
        end
    elseif Rnd == Rnum
        for (Inum, Idenom) in zip(Rnum, Rdenom)
            @inbounds numdenom[Inum] = NumDenom(num[Inum], denom[Idenom])
        end
    else
        for (Ind, Inum, Idenom) in zip(Rnd, Rnum, Rdenom)
            @inbounds numdenom[Ind] = NumDenom(num[Inum], denom[Idenom])
        end
    end
    numdenom
end

function _packnd!(numdenom::CenterIndexedArray, num::CenterIndexedArray, denom::CenterIndexedArray)
    @simd for I in eachindex(num)
        @inbounds numdenom[I] = NumDenom(num[I], denom[I])
    end
    numdenom
end

# The next are mostly used just for testing
"""
`mms = mismatcharrays(nums, denoms)` packs array-of-arrays num/denom pairs as an array-of-MismatchArrays.

`mms = mismatcharrays(nums, denom)`, for `denom` a single array, uses the same `denom` array for all `nums`.
"""
function mismatcharrays(nums::AbstractArray{A}, denom::AbstractArray{T}) where {A<:AbstractArray,T<:Number}
    first = true
    local mms
    for i in eachindex(nums)
        num = nums[i]
        mm = MismatchArray(num, denom)
        if first
            mms = Array{typeof(mm)}(undef, size(nums))
            first = false
        end
        mms[i] = mm
    end
    mms
end

function mismatcharrays(nums::AbstractArray{A1}, denoms::AbstractArray{A2}) where {A1<:AbstractArray,A2<:AbstractArray}
    size(nums) == size(denoms) || throw(DimensionMismatch("nums and denoms arrays must have the same number of apertures"))
    first = true
    local mms
    for i in eachindex(nums, denoms)
        mm = MismatchArray(nums[i], denoms[i])
        if first
            mms = Array{typeof(mm)}(undef, size(nums))
            first = false
        end
        mms[i] = mm
    end
    mms
end

"""
`num, denom = separate(mm)` splits an `AbstractArray{NumDenom}` into separate
numerator and denominator arrays.
"""
function separate(data::AbstractArray{NumDenom{T}}) where T
    num = Array{T}(undef, size(data))
    denom = similar(num)
    for I in eachindex(data)
        nd = data[I]
        num[I] = nd.num
        denom[I] = nd.denom
    end
    num, denom
end

function separate(mm::MismatchArray)
    num, denom = separate(mm.data)
    CenterIndexedArray(num), CenterIndexedArray(denom)
end

function separate(mma::AbstractArray{M}) where M<:MismatchArray
    T = eltype(eltype(M))
    nums = Array{CenterIndexedArray{T,ndims(M)}}(undef, size(mma))
    denoms = similar(nums)
    for (i,mm) in enumerate(mma)
        nums[i], denoms[i] = separate(mm)
    end
    nums, denoms
end

"""
    r = ratio(nd::NumDenom, thresh, fillval=NaN)

Return `nd.num/nd.denom`, unless `nd.denom < thresh`, in which case return `fillval` converted
to the same type as the ratio.
Choosing a `thresh` of zero will always return the ratio.
"""
@inline function ratio(nd::NumDenom{T}, thresh, fillval=convert(T,NaN)) where {T}
    r = nd.num/nd.denom
    return nd.denom < thresh ? oftype(r, fillval) : r
end
ratio(r::Real, thresh, fillval=NaN) = r

(::Type{M})(::Type{T}, dims::Dims) where {M<:MismatchArray,T} = CenterIndexedArray{NumDenom{T}}(undef, dims)
(::Type{M})(::Type{T}, dims::Integer...) where {M<:MismatchArray,T} = CenterIndexedArray{NumDenom{T}}(undef, dims)

function Base.copyto!(M::MismatchArray, nd::Tuple{AbstractArray, AbstractArray})
    num, denom = nd
    size(M) == size(num) == size(denom) || error("all sizes must match")
    for (IM, Ind) in zip(eachindex(M), eachindex(num))
        M[IM] = NumDenom(num[Ind], denom[Ind])
    end
    M
end


#### Utility functions ####

"""
`index = indmin_mismatch(numdenom, thresh)` returns the location of
the minimum value of what is effectively `num./denom`.  However, it
considers only those points for which `denom .> thresh`; moreover, it
will never choose an edge point.  `index` is a CartesianIndex into the
arrays.
"""
function indmin_mismatch(numdenom::MismatchArray{NumDenom{T},N}, thresh::Real) where {T,N}
    imin = CartesianIndex(ntuple(d->0, Val(N)))
    rmin = typemax(T)
    threshT = convert(T, thresh)
    @inbounds for I in CartesianIndices(map(trimedges, axes(numdenom)))
        nd = numdenom[I]
        if nd.denom > threshT
            r = nd.num/nd.denom
            if r < rmin
                imin = I
                rmin = r
            end
        end
    end
    return imin
end

function indmin_mismatch(r::CenterIndexedArray{T,N}) where {T<:Number,N}
    imin = CartesianIndex(ntuple(d->0, Val(N)))
    rmin = typemax(T)
    @inbounds for I in CartesianIndices(map(trimedges, axes(r)))
        rval = r[I]
        if rval < rmin
            imin = I
            rmin = rval
        end
    end
    return imin
end

trimedges(r::AbstractUnitRange) = first(r)+1:last(r)-1

### Miscellaneous

"""
`datahp = highpass([T], data, sigma)` returns a highpass-filtered
version of `data`, with all negative values truncated at 0.  The
highpass is computed by subtracting a lowpass-filtered version of
data, using Gaussian filtering of width `sigma`.  As it is based on
`Image.jl`'s Gaussian filter, it gracefully handles `NaN` values.

If you do not wish to highpass-filter along a particular axis, put
`Inf` into the corresponding slot in `sigma`.

You may optionally specify the element type of the result, which for
`Integer` or `FixedPoint` inputs defaults to `Float32`.
"""
function highpass(::Type{T}, data::AbstractArray, sigma) where T
    if any(isinf, sigma)
        datahp = convert(Array{T,ndims(data)}, data)
    else
        datahp = data - imfilter(T, data, KernelFactors.IIRGaussian(T, (sigma...,)), NA())
    end
    datahp[datahp .< 0] .= 0  # truncate anything below 0
    datahp
end
highpass(data::AbstractArray{T}, sigma) where {T<:AbstractFloat} = highpass(T, data, sigma)
highpass(data::AbstractArray, sigma) = highpass(Float32, data, sigma)

"""
`pp = PreprocessSNF(bias, sigmalp, sigmahp)` constructs an object that
can be used to pre-process an image as `pp(img)`. The "SNF" part of
the name means "shot-noise filtered," meaning that this preprocessor
is specifically designed for situations in which you are dominated by
shot noise (i.e., from photon-counting statistics).

The processing is of the form
```
    imgout = bandpass(√max(0,img-bias))
```
i.e., the image is bias-subtracted, square-root transformed (to turn
shot noise into constant variance), and then band-pass filtered using
Gaussian filters of width `sigmalp` (for the low-pass) and `sigmahp`
(for the high-pass).  You can pass `sigmalp=zeros(n)` to skip low-pass
filtering, and `sigmahp=fill(Inf, n)` to skip high-pass filtering.
"""
mutable struct PreprocessSNF  # Shot-noise filtered
    bias::Float32
    sigmalp::Vector{Float32}
    sigmahp::Vector{Float32}
end
# PreprocessSNF(bias::T, sigmalp, sigmahp) = PreprocessSNF{T}(bias, T[sigmalp...], T[sigmahp...])

function preprocess(pp::PreprocessSNF, A::AbstractArray)
    Af = sqrt_subtract_bias(A, pp.bias)
    imfilter(highpass(Af, pp.sigmahp), KernelFactors.IIRGaussian((pp.sigmalp...,)), NA())
end
(pp::PreprocessSNF)(A::AbstractArray) = preprocess(pp, A)
# For SubArrays, extend to the parent along any non-sliced
# dimension. That way, we keep any information from padding.
function (pp::PreprocessSNF)(A::SubArray)
    Bpad = preprocess(pp, paddedview(A))
    trimmedview(Bpad, A)
end
# ImageMeta method is defined under @require in __init__

function sqrt_subtract_bias(A, bias)
#    T = typeof(sqrt(one(promote_type(eltype(A), typeof(bias)))))
    T = Float32
    out = Array{T}(undef, size(A))
    for I in eachindex(A)
        @inbounds out[I] = sqrt(max(zero(T), convert(T, A[I]) - bias))
    end
    out
end


"""
`Apad = paddedview(A)`, for a SubArray `A`, returns a SubArray that
extends to the full parent along any non-sliced dimensions of the
parent. See also [`trimmedview`](@ref).
"""
paddedview(A::SubArray) = _paddedview(A, (), (), A.indices...)
_paddedview(A::SubArray{T,N,P,I}, newindexes, newsize) where {T,N,P,I} =
    SubArray(A.parent, newindexes)
@inline function _paddedview(A, newindexes, newsize, index, indexes...)
    d = length(newindexes)+1
    _paddedview(A, (newindexes..., pdindex(A.parent, d, index)), pdsize(A.parent, newsize, d, index), indexes...)
end
pdindex(A, d, i::Base.Slice) = i
pdindex(A, d, i::Real) = i
pdindex(A, d, i::UnitRange) = 1:size(A,d)
pdindex(A, d, i) = error("Cannot pad with an index of type ", typeof(i))

pdsize(A, newsize, d, i::Base.Slice) = tuple(newsize..., size(A,d))
pdsize(A, newsize, d, i::Real) = newsize
pdsize(A, newsize, d, i::UnitRange) = tuple(newsize..., size(A,d))

"""
`B = trimmedview(Bpad, A::SubArray)` returns a SubArray `B` with
`axes(B) = axes(A)`. `Bpad` must have the same size as `paddedview(A)`.
"""
function trimmedview(Bpad, A::SubArray)
    ndims(Bpad) == ndims(A) || throw(DimensionMismatch("dimensions $(ndims(Bpad)) and $(ndims(A)) of Bpad and A must match"))
    _trimmedview(Bpad, A.parent, 1, (), A.indices...)
end
_trimmedview(Bpad, P, d, newindexes) = view(Bpad, newindexes...)
@inline _trimmedview(Bpad, P, d, newindexes, index::Real, indexes...) =
    _trimmedview(Bpad, P, d+1, newindexes, indexes...)
@inline function _trimmedview(Bpad, P, d, newindexes, index, indexes...)
    dB = length(newindexes)+1
    Bsz = size(Bpad, dB)
    Psz = size(P, d)
    Bsz == Psz || throw(DimensionMismatch("dimension $dB of Bpad has size $Bsz, should have size $Psz"))
    _trimmedview(Bpad, P, d+1, (newindexes..., index), indexes...)
end

# For faster and type-stable slicing
struct ColonFun end
ColonFun(::Int) = Colon()

function __init__()
    @require ImageMetadata="bc367c6b-8a6b-528e-b4bd-a4b897500b49" begin
        (pp::PreprocessSNF)(A::ImageMetadata.ImageMeta) = ImageMetadata.shareproperties(A, pp(ImageMetadata.arraydata(A)))
    end
end

include("deprecations.jl")

end
