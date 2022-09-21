# resample pixels on floe boundary for psi-s curve


"""
    resample_boundary(bd_points::Vector{CartesianIndex}, reduc_factor::Int64=2, bd::String="natural")

Get a uniform set of resampled boundary points from `bd_points` using cubic splines with specified boundary conditions

The resampled set of points is obtained using parametric interpolation of the points in `bd_points`. It is assumed that the separation between a pair of adjacent points is 1.

# Arguments
- `bd_points`: Sequetial set of boundary points for the object of interest
- `reduc_factor`: Factor by which to reduce the number of points in `bd_points` (2 by default)
-`bd`: Boundary condition, either 'natural' (default) or 'periodic'

See also [`bwtraceboundary`](@ref)

# Example

```jldoctest; setup = :(using IceFloeTracker)
julia> A = zeros(Int, 13, 16); A[2:6, 2:6] .= 1; A[4:8, 7:10] .= 1; A[10:12,13:15] .= 1; A[10:12,3:6] .= 1;
julia> A
13×16 Matrix{Int64}:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0
 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0
 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0
 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

julia> boundary = bwtraceboundary(A);

julia> boundary[3]
9-element Vector{CartesianIndex}:
 CartesianIndex(10, 13)
 CartesianIndex(11, 13)
 CartesianIndex(12, 13)
 CartesianIndex(12, 14)
 CartesianIndex(12, 15)
 CartesianIndex(11, 15)
 CartesianIndex(10, 15)
 CartesianIndex(10, 14)
 CartesianIndex(10, 13)

 julia> resample_boundary(boundary[3])
4×2 Matrix{Float64}:
 10.0     13.0
 12.0357  13.5859
 10.5859  15.0357
 10.0     13.0
"""
function resample_boundary(bd_points::Vector{<:CartesianIndex}, reduc_factor::Int64=2, bd::String="natural")
    # check boundary conditions
    if bd == "natural"
        BD = Natural(OnGrid())
    elseif bd == "periodic"
        BD = Periodic(OnGrid())
    end
    
    # reparemetrize using arclength
    s_in = range(0,1,length(bd_points))

    # arclengths to resample
    s_out = range(0,1,length(bd_points) ÷ reduc_factor)

    # collect data in bd_points for interpolant
    A = getindex.(bd_points, [1 2])
    
    # build interpolant generator
    itp = scale(interpolate(A, (BSpline(Cubic(BD)), NoInterp())), s_in, 1:2)

    # get resampled data
    xs, ys = [itp(s,1) for s in s_out], [itp(s,2) for s in s_out]
    
    return [xs ys]
end
