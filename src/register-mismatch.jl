# Pkg.add("RegisterMismatch"); Pkg.add("RegisterQD")

"""
    mismatch(fixed::AbstractArray,
                  moving::AbstractArray,
                  mxshift::Tuple{Int64, Int64}=(10,10),
                  mxrot::Float64=pi/4;
                  kwargs...
                  )                   

Estimate a rigid transformation (translation + rotation) that minimizes the 'mismatch' of aligning `moving` with `fixed` using the [QuadDIRECT algorithm](https://github.com/timholy/QuadDIRECT.jl#readme).

Returns a pair with the mismatch score `mm` and the associated transformation `tfm`.

# Arguments
- `fixed`,`moving`: images to align via a rigid transformation
- `mxshift`: maximum allowed translation in units of array indices (default set to `(10,10)`)
- `mxrot`: maximum allowed rotation in radians (default set to `π/4`)
- `thresh`: minimum sum-of-squared-intensity overlap between the images (default is 10% of the sum-of-squared-intensity of `fixed`)
- `kwargs`: other arguments such as `tol`, `ftol`, and `fvalue` (see [QuadDIRECT.analyze](https://github.com/timholy/QuadDIRECT.jl/blob/d6170f14a49f57552c59c9b4533a4b75a3ab3c45/src/algorithm.jl#L459) for details)

```
"""
function mismatch(
    fixed::AbstractArray,
    moving::AbstractArray,
    mxshift::Tuple{Int64,Int64}=(10, 10),
    mxrot::Float64=pi / 4;
    kwargs...,
)
    tfm, mm = qd_rigid(
        centered(fixed), centered(moving), mxshift, mxrot; print_interval=typemax(Int)
    )

    return mm, tfm
end
