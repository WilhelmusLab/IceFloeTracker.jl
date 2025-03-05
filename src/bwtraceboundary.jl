# Adapted from https://juliaimages.org/latest/examples/contours/contour_detection/#Contour-Detection-and-Drawing

# Boundary tracing algorithm from
# Suzuki, Satoshi. "Topological structural analysis of digitized binary images by border following." Computer vision, graphics, and image processing 30.1 (1985): 32-46.

"""
    bwtraceboundary(image::Union{Matrix{Int64},Matrix{Float64},T};
                    P0::Union{Tuple{Int,Int},CartesianIndex{2},Nothing}=nothing,
                    closed::Bool=true) where T<:AbstractMatrix{Bool}

Trace the boundary of objects in `image`

Background pixels are represented as zero. The algorithm traces the boundary _counterclockwise and an initial point `P0` can be specified. If more than one boundary is detected and an initial point is provided, the boundary that contains this point is returned as a vector of CartesianIndex types. Otherwise an array of vectors is returned with all the detected boundaries in `image`.

# Arguments
- `image`: image, preferably binary with one single object, whose objects' boundaries are to be traced.
- `P0`: initial point of a target boundary.
- `closed`: if `true` (default) makes the inital point of a boundary equal to the last point.

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

julia> boundary = IceFloeTracker.bwtraceboundary(A);

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
```
"""
function bwtraceboundary(
    image::Union{Matrix{UInt8},Matrix{Int64},Matrix{Float64},T};
    P0::Union{Tuple{Int,Int},CartesianIndex{2},Nothing}=nothing,
    closed::Bool=true,
) where {T<:AbstractArray{Bool,2}}
    if typeof(image[1]) != Float64
        image = Float64.(image)
    end

    nbd = 1

    contour_list = Vector{typeof(CartesianIndex[])}()
    done = falses(8)

    # Clockwise Moore neighborhood.
    dir_delta = [
        CartesianIndex(-1, 0),
        CartesianIndex(-1, 1),
        CartesianIndex(0, 1),
        CartesianIndex(1, 1),
        CartesianIndex(1, 0),
        CartesianIndex(1, -1),
        CartesianIndex(0, -1),
        CartesianIndex(-1, -1),
    ]

    height, width = size(image)

    for i in 1:height
        for j in 1:width
            is_outer = (image[i, j] == 1 && (j == 1 || image[i, j - 1] == 0)) ## 1 (a)
            is_hole = (image[i, j] >= 1 && (j == width || image[i, j + 1] == 0))

            if is_outer || is_hole
                # 2
                border = CartesianIndex[]

                from = CartesianIndex(i, j)

                if is_outer
                    nbd += 1
                    from -= CartesianIndex(0, 1)
                else
                    nbd += 1
                    from += CartesianIndex(0, 1)
                end

                p0 = CartesianIndex(i, j)
                _detect_move!(image, p0, from, nbd, border, done, dir_delta) ## 3

                if isempty(border)
                    push!(border, p0)
                    image[p0] = -nbd
                end
                push!(contour_list, border)
            end
        end
    end

    # make contour start at provided P0

    if .!isnothing(P0)
        P0 = CartesianIndex(P0)

        # check P0 is in a contour
        test_p0, idx = isincountourlist(P0, contour_list)

        # if so, make the contour start at P0
        if test_p0
            target_contour = contour_list[idx]
            target_contour = makecontourstartatP(P0, target_contour)
            if closed
                return push!(target_contour, target_contour[1])
            else
                return target_contour
            end
        else
            @warn "Point at $(Tuple(P0)) not found in any contour. Returning all found countours."
        end
    end

    # if P0 not in any contour but want each sequence closed
    if closed
        for seq in contour_list
            push!(seq, seq[1])
        end
    end

    # unpack contour_list in case only one contour is found
    return contour_list
end

"""
Make contour start at point P by permuting the elements in contour.
"""
function makecontourstartatP(P::CartesianIndex{2}, contour::Vector{CartesianIndex})
    iout = 1 # initialize variable for use in for loop

    # find index of P in contour
    for (i, cidx) in enumerate(contour)
        if P == cidx
            iout = i
            break
        end
    end

    if iout == 1
        return contour
    else
        # return correct permutation of contour
        return [contour[iout:end]; contour[1:(iout - 1)]]
    end
end

"""
Check `P` is in countour list if so return the index of the contour that contains `P`, otherwise return false.
"""
function isincountourlist(
    P::Union{CartesianIndex{2},Tuple{Int64,Int64}},
    contour_list::Vector{Vector{CartesianIndex}},
)
    if typeof(P) <: Tuple
        P = CartesianIndex(P)
    end

    for (i, contour) in enumerate(contour_list)
        if P in contour
            return (true, i)
        end
    end

    return (false, 0)
end

# """
# Get index in Moore neigborhood representing the direction from the `from` pixel coords to the `to` pixel coords (see definition of dir_delta below).

# ## Clockwise Moore neighborhood.
# dir_delta = [CartesianIndex(-1, 0), CartesianIndex(-1, 1), CartesianIndex(0, 1), CartesianIndex(1, 1), CartesianIndex(1, 0), CartesianIndex(1, -1), CartesianIndex(0, -1), CartesianIndex(-1,-1)]

# """
function _from_to(
    from::CartesianIndex, to::CartesianIndex, dir_delta::Vector{CartesianIndex{2}}
)
    delta = to - from
    return findfirst(x -> x == delta, dir_delta)
end

# """
# Workhorse function: Get all pixel coords for detected border.
# """
function _detect_move!(
    image::Matrix{Float64},
    p0::CartesianIndex{2},
    p2::CartesianIndex{2},
    nbd::Int,
    border::Vector{CartesianIndex},
    done::BitVector,
    dir_delta::Vector{CartesianIndex{2}},
)
    dir = _from_to(p0, p2, dir_delta)
    moved = _clockwise(dir)
    p1 = CartesianIndex(0, 0)

    while moved != dir ## 3.1
        newp = move(p0, image, moved, dir_delta)
        if newp[1] != 0
            p1 = newp
            break
        end
        moved = _clockwise(moved)
    end

    if p1 == CartesianIndex(0, 0)
        return nothing
    end

    p2 = p1 ## 3.2
    p3 = p0 ## 3.2
    done .= false

    while true
        dir = _from_to(p3, p2, dir_delta)
        moved = _counterclockwise(dir)
        p4 = CartesianIndex(0, 0) # initialize p4
        done .= false

        while p4[1] == 0 ## 3.3: "Examine N(p3) for a nonzero pixel, name first nonzero pixel as p4"
            p4 = move(p3, image, moved, dir_delta)
            done[moved] = true
            moved = _counterclockwise(moved)
        end

        push!(border, p3) ## 3.4

        if p3[1] == size(image, 1) || done[3]
            image[p3] = -nbd
        elseif image[p3] == 1
            image[p3] = nbd
        end

        if (p4 == p0 && p3 == p1) ## 3.5
            break
        end
        p2 = p3
        p3 = p4
    end
end

# """
# Make a clockwise turn from the `dir` direction
# """
function _clockwise(dir::Int)
    return (dir) % 8 + 1
end

# """
# Make a counterclockwise turn from the `dir` direction
# """
function _counterclockwise(dir::Int)
    return (dir + 6) % 8 + 1
end

"""
move from current pixel to the next in given direction
"""
function move(
    pixel::CartesianIndex{2},
    image::Matrix{Float64},
    dir::Int,
    dir_delta::Vector{CartesianIndex{2}},
)
    newp = pixel + dir_delta[dir]
    height, width = size(image)
    if (0 < newp[1] <= height) && (0 < newp[2] <= width) # check move is within image domain
        if image[newp] != 0 # and that it is not a background point
            return newp
        end
    end
    return CartesianIndex(0, 0)
end
