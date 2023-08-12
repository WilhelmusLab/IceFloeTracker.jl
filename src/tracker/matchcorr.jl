"""
    matchcorr(
    f1::T,
    f2::T,
    Δt::F,
    mxrot::S=10,
    psi::F=0.95,
    sz::S=16,
    comp::F=0.25,
    mm::F=0.22
    )
    where {T<:AbstractArray{Bool,2},S<:Int64,F<:Float64}

Compute the mismatch `mm` and psi-s-correlation `c` for floes with masks `f1` and `f2`.

The criteria for floes to be considered equivalent is as follows:
    - `c` greater than `mm` 
    - `_mm` is less than `mm`

A pair of `NaN` is returned for cases for which one of their mask dimension is too small or their sizes are not comparable.

# Arguments
- `f1`: mask of floe 1
- `f2`: mask of floe 2
- `Δt`: time difference between floes
- `mxrot`: maximum rotation (in degrees) allowed between floesn (default: 10)
- `psi`: psi-s-correlation threshold (default: 0.95)
- `sz`: size threshold (default: 16)
- `comp`: size comparability threshold (default: 0.25)
- `mm`: mismatch threshold (default: 0.22)
"""
function matchcorr(
    f1::T, f2::T, Δt::F; mxrot::S=10, psi::F=0.95, sz::S=16, comp::F=0.25, mm::F=0.22
) where {T<:AbstractArray{Bool,2},S<:Int64,F<:Float64}

    # check if the floes are too small and size are comparable
    _sz = size.([f1, f2])
    if (any([(_sz...)...] .< sz) || getsizecomparability(_sz...) > comp)
        return (mm=NaN, c=NaN)
    end

    _psi = buildψs.([f1, f2])
    c = corr(_psi...)

    if c < psi
        @warn "correlation too low, c: $c"
        return (mm=NaN, c=NaN)
    else
        return (mm=0.0, c=c)
    end

    # check if the time difference is too large or the rotation is too large
    _mm, rot = mismatch(f1, f2; mxrot=deg2rad(mxrot))
    if all([Δt < 300, rot > mxrot]) || _mm > mm
        @warn "time difference too small for a large rotation or mismatch too large\nmm: $mm, rot: $rot"
        return (mm=NaN, c=NaN)
    end
    if mm < 0.1
        mm = 0.0
    end
    return (mm=mm, c=c)
end

"""
    getsizecomparability(s1, s2)

Check if the size of two floes `s1` and `s2` are comparable. The size is defined as the product of the floe dimensions.

# Arguments
- `s1`: size of floe 1
- `s2`: size of floe 2
"""
function getsizecomparability(s1::T, s2::T) where {T<:Tuple{Int64,Int64}}
    a1 = *(s1...)
    a2 = *(s2...)
    return abs(a1 - a2) / a1
end

"""
    corr(f1,f2)

Return the normalized cross-correlation between the psi-s curves `p1` and `p2`.
"""
function corr(p1::T, p2::T) where {T<:AbstractArray}
    cc, _ = maximum.(IceFloeTracker.crosscorr(p1, p2; normalize=true))
    return cc
end
