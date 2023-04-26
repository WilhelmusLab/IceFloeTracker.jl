"""
    matchcorr(
    f1, f2, Δt, mxrot=10, psi_s_thresh=0.95, sz_thresh=16, comp_tresh=0.25, mm_thresh=0.22
)


    matchcorr(f1, f2, Δt, mxrot, psi_s_thresh, sz_thresh, comp_tresh)

Compute the mismatch `mm` and psi-correlation `c` for floes with masks `f1` and `f2`.

The criteria for floes to be considered equivalent is as follows:
    - `c` greater than `mm_thresh` 
    - `mm` is less than `mm_thresh`

A pair of `NaN` is returned for cases for which one of their mask dimension is too small or their sizes are not comparable.

# Arguments
- `f1`: mask of floe 1
- `f2`: mask of floe 2
- `Δt`: time difference between floes
- `mxrot`: maximum rotation allowed between floes
- `psi_thresh`: psi-correlation threshold
- `sz_thresh`: size threshold
- `comp_tresh`: size comparability threshold
- `mm_thresh`: mismatch threshold
"""
function matchcorr(
    f1, f2, Δt, mxrot=10, psi_thresh=0.95, sz_thresh=16, comp_tresh=0.25, mm_thresh=0.22
)

    # check if the floes are too small and size are comparable
    sz = size.([f1, f2])
    (any([(sz...)...] .< sz_thresh) || getsizecomparability(sz...) > comp_tresh) &&
        return (mm=NaN, c=NaN)

    psi = buildψs.([f1, f2])
    c = corr(psi...)

    c < psi_thresh && return (mm=NaN, c=NaN)

    # check if the time difference is too large or the rotation is too large
    mm, rot = mismatch(f1, f2; mxrot=deg2rad(mxrot))
    any([Δt < 300, rot > mxrot, mm > mm_thresh]) && return (mm=NaN, c=NaN)

    mm < 0.1 && return (mm=0, c=c)

    return (mm=mm, c=c)
end

"""
    getsizecomparability(s1, s2)

Check if the size of two floes `s1` and `s2` are comparable. The size is defined as the product of the floe dimensions.

# Arguments
- `s1`: size of floe 1
- `s2`: size of floe 2
"""
function getsizecomparability(s1, s2)
    a1 = *(s1...)
    a2 = *(s2...)
    return abs(a1 - a2) / a1
end

"""
    corr(f1,f2)

Return the correlation between the psi-s curves `p1` and `p2`.
"""
function corr(p1, p2)
    cc, _ = maximum.(IceFloeTracker.crosscorr(p1, p2; normalize=true))
    return cc
end
