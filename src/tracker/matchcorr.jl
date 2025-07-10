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
- `Δt`: time difference between floes in minutes
- `mxrot`: maximum rotation (in degrees) allowed between floes (default: 10)
- `psi`: psi-s-correlation threshold (default: 0.95)
- `sz`: size threshold (default: 16)
- `comp`: size comparability threshold (default: 0.25)
- `mm`: mismatch threshold (default: 0.22)
"""
function matchcorr(
    f1::T, f2::T, Δt::F; mxrot::S=10, psi::F=0.95, sz::S=16, comp::F=0.25, mm::F=0.22
) where {T<:AbstractArray{Bool,2},S<:Int64,F<:Float64}

    # tbd: add to function signature 
    # setting to 10 days for now
    max_dt_minutes = 10*24*60

    # tbd: add to function signature
    # confidence level critical number: default is for 95%
    cn = 1.96
    
    # check if the floes are too small and size are comparable
    _sz = size.([f1, f2])
    if (any([(_sz...)...] .< sz) || getsizecomparability(_sz...) > comp)
        return (mm=NaN, c=NaN)
    end

    _psi = buildψs.([f1, f2])
    r = round(corr(_psi...), digits=3)
    
    # confidence interval for Pearson correlation coefficient
    z = 0.5*log((1 + r)/(1 - r))
    n = minimum(length.(_psi))
    sigma_z = sqrt(1/(n - 3))
    zlow = z - cn*sigma_z
    zhigh = z + cn*sigma_z
    rlow = (zlow - 1)/(zlow + 1)
    rhigh = (zhigh - 1)/(zhigh + 1)
    
    if r < psi
        @warn "correlation too low, r: $r"
        return (mm=NaN, c=r)
    # dmw: am I wrong in thinking this avoids calculating mismatch completely?
    # else
    #     return (mm=0.0, c=c)
    end

    # check if the time difference is too large or the rotation is too large
    _mm, rot = mismatch(f1, f2; mxrot=deg2rad(mxrot))
    if all([Δt < max_dt_minutes, rot > mxrot]) || _mm > mm
        @warn "time difference too small for a large rotation or mismatch too large\nmm: $mm, rot: $rot"
        return (mm=NaN, c=NaN)
    end

    # dmw: why would we erase the mismatch measurement?
    # if mm < 0.1
    #     mm = 0.0
    # end
    
    return (mm=mm, c=r, ci=(rlow, rhigh))
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
# dmw: Can we rename this to avoid confusion with standard Pearson correlation?
function corr(p1::T, p2::T) where {T<:AbstractArray}
    cc, _ = maximum.(IceFloeTracker.crosscorr(p1, p2; normalize=true))
    return cc
end
