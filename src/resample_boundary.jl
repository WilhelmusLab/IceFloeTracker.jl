# resample pixels on floe boundary for psi-s curve

"""
    resample_boundary(bd_points::Vector{CartesianIndex}, reduc_factor::Int64=2)

Get a uniform set of resampled boundary points from `bd_points` using cubic splines with 'natural' boundary conditions

The resampled set of points is obtained using parametric interpolation of the points in `bd_points`. It is assumed that the separation between a pair of adjacent points is 1.

# Arguments
- `bd_points`: Sequetial set of boundary points for the object of interest
- `reduc_factor`: factor by which to reduce the number of points in `bd_points`

See also [`bwtraceboundary`](@ref)
"""
function resample_boundary(bd_points::Vector{CartesianIndex}, reduc_factor::Int64=2)
    
    # reparemetrize using arclength
    s_in = range(0,1,length(bd_points))

    # arclengths to resample
    s_out = range(0,1,length(bd) ÷ reduc_factor)

    # collect data in bd_points for interpolant
    A = getindex.(bd, [1 2])
    
    # build interpolant generator
    itp = scale(interpolate(A, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), s_in, 1:2)

    # get resampled data
    xs, ys = [itp(t,1) for t in s_out], [itp(t,2) for t in s_out]
    
    return (xs,ys)
end
