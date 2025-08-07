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

Compute the mismatch `sd` and psi-s-correlation `c` for floes with masks `f1` and `f2`.

The criteria for floes to be considered good matches is as follows:
    - `c` greater than `psi` 
    - sd is less than `mm`

A pair of `NaN` is returned for cases for which one of their mask dimension is too small or their sizes are not comparable.

# Arguments
- `f1`: mask of floe 1
- `f2`: mask of floe 2
- `Δt`: time difference between floes in minutes
- `rot_stepsize`: step size in degrees for registration
- `mxrot`: maximum rotation (in degrees) allowed between floes (default: 10)
- `psi`: psi-s-correlation threshold (default: 0.95)
- `sz`: size threshold (default: 16)
- `comp`: size comparability threshold (default: 0.25)
- `mm`: mismatch threshold (default: 0.22)
"""
function matchcorr(
    f1::T, f2::T, Δt::F; rot_stepsize::S=1, mxrot::S=10, psi::F=0.95, sz::S=16, comp::F=0.25, mm::F=0.22 # update varnames
) where {T<:AbstractArray{Bool,2},S<:Int64,F<:Float64}

    # tbd: add to function signature 
    # setting to 10 days for now
    max_dt_minutes = 10*24*60

    # tbd: add to function signature
    # confidence level critical number: default is for 95%
    cn = 1.96
    
    # check if the floes are too small and size are comparable
    # dmw: this step is redundant, since we are checking the absolute difference ratios earlier in the tracker
    _sz = size.([f1, f2])
    if (any([(_sz...)...] .< sz) || getsizecomparability(_sz...) > comp)
        return (shape_difference=NaN,
            psi_s_correlation=NaN,
            rotation=NaN,
            corr_ci=(NaN, NaN),
            sd_ci=(NaN, NaN),
            rotation_ci=(NaN, NaN))
    end

    # dmw: ψ-s curves are also in the region props table, we could use those instead of re-calculating
    _psi = buildψs.([f1, f2])
    r = round(normalized_maximum_crosscorrelation(_psi...), digits=3)
    
    # confidence interval for Pearson correlation coefficient
    z = 0.5*log((1 + r)/(1 - r))
    n = minimum(length.(_psi))
    sigma_z = sqrt(1/(n - 3))
    zlow = z - cn*sigma_z
    zhigh = z + cn*sigma_z
    rlow = round((exp(2*zlow) - 1)/(exp(2*zlow) + 1), digits=3)
    rhigh = round((exp(2*zhigh) - 1)/(exp(2*zhigh) + 1), digits=3)
    
    if r < psi
        @warn "correlation too low, r: $r"
        return (shape_difference=NaN,
            psi_s_correlation=r,
            rotation=NaN,
            corr_ci=(rlow, rhigh),
            sd_ci=(NaN, NaN),
            rotation_ci=(NaN, NaN))
    end

    # warn if the time difference is too large or the rotation is too large
    # dmw: think through this section, should the matchcorr function be constraining these
    # terms or should the constraints occur in the tracker?
    # dmw: it looks like the mismatch function isn't registering the stepsize parameter.
    _mm, rot = mismatch(f1, f2, mxrot, rot_stepsize)
    _mm = round(_mm, digits=3)
    rot = round(rot, digits=3)
    if all([Δt < max_dt_minutes, rot > mxrot]) || _mm > mm
        @warn "time difference too small for a large rotation or mismatch too large\nmm: $_mm (threshold $mm), rot: $rot (threshold $mxrot)"
    end

    return (shape_difference=_mm,
            psi_s_correlation=r,
            rotation=rot,
            corr_ci=(rlow, rhigh),
            sd_ci=(NaN, NaN),
            rotation_ci=(NaN, NaN))
end

"""
    getsizecomparability(s1, s2)

Check if the size of the bounding boxes for two floes `s1` and `s2` are comparable. The comparibility is based on the absolute difference ratio.

# Arguments
- `s1`: size of floe 1
- `s2`: size of floe 2
"""
function getsizecomparability(s1::T, s2::T) where {T<:Tuple{Int64,Int64}}
    a1 = *(s1...)
    a2 = *(s2...)
    return abs(a1 - a2) / (0.5 * (a1 + a2))
end

"""
    corr(f1,f2)

Return the normalized cross-correlation between the psi-s curves `p1` and `p2`.
"""
# dmw: Renamed this to avoid confusion with standard Pearson correlation
# to do: check if the curves need to be psi-s curves, or if any vector is ok
function normalized_maximum_crosscorrelation(p1::T, p2::T) where {T<:AbstractArray}
    cc, _ = maximum.(IceFloeTracker.crosscorr(p1, p2; normalize=true))
    return cc
end