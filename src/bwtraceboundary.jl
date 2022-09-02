# Adapted from https://juliaimages.org/latest/examples/contours/contour_detection/#Contour-Detection-and-Drawing

# Boundary tracing algorithm from
# Suzuki, Satoshi. "Topological structural analysis of digitized binary images by border following." Computer vision, graphics, and image processing 30.1 (1985): 32-46.


# To do: create auxiliary methods

"""
    bwtraceboundary(image::AbstractArray{2},
                    P0::Union{Tuple{Int,Int},CartesianIndex{2}}=nothing,
                    closed::Bool=true)

Trace the boundary of objects in `image` 

Background pixels are represented as zeros. The algorithm traces the boundary counterclockwise and an initial point `P0` can be specified. If more than one boundary is detected and an initial point is provided, the boundary that contains this point is returned.

"""
function bwtraceboundary(image,
                         P0::Union{Tuple{Int,Int},CartesianIndex{2},Any}=nothing,
                         closed::Union{Bool,Any}=true)
    nbd = 1
    lnbd = 1
    image = Float64.(image)
    contour_list =  Vector{typeof(CartesianIndex[])}()
    done = [false, false, false, false, false, false, false, false]

    # Clockwise Moore neighborhood.
    dir_delta = [CartesianIndex(-1, 0) , CartesianIndex(-1, 1), CartesianIndex(0, 1), CartesianIndex(1, 1), CartesianIndex(1, 0), CartesianIndex(1, -1), CartesianIndex(0, -1), CartesianIndex(-1,-1)]

    height, width = size(image)

    for i=1:height
        lnbd = 1
        for j=1:width
            fji = image[i, j]
            is_outer = (image[i, j] == 1 && (j == 1 || image[i, j-1] == 0)) ## 1 (a)
            is_hole = (image[i, j] >= 1 && (j == width || image[i, j+1] == 0))

            if is_outer || is_hole
                # 2
                border = CartesianIndex[]

                from = CartesianIndex(i, j)

                if is_outer
                    nbd += 1
                    from -= CartesianIndex(0, 1)

                else
                    nbd += 1
                    if fji > 1
                        lnbd = fji
                    end
                    from += CartesianIndex(0, 1)
                end

                p0 = CartesianIndex(i,j)
                detect_move(image, p0, from, nbd, border, done, dir_delta) ## 3
                if isempty(border) ##TODO
                    push!(border, p0)
                    image[p0] = -nbd
                end
                push!(contour_list, border)
            end
            if fji != 0 && fji != 1
                lnbd = abs(fji)
            end

        end
    end

    # make contour start at provided P0
    
    if .!isnothing(P0)

        P0 = CartesianIndex(P0)

        # check P0 is in a contour
        test_p0 = isincountourlist(P0,contour_list)


        # if so, make the contour start at P0
        if test_p0 != false
            target_contour = contour_list[test_p0]
            target_contour = makecontourstartatP(P0, target_contour)
            if closed
                return push!(target_contour,target_contour[1])
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
    
    return contour_list
end

"""
Make contour start at point P by permuting the elements in contour.
"""
function makecontourstartatP(P, contour)

    iout = 1
    # contourout = nothing

    # find index of P in contour
    for (i, cidx) in enumerate(contour)
        if P == cidx
            iout = i
            break
        end
    end

    if iout == 1
        return contour
    elseif iout > 1
        # return correct permutation of contour
        return [contour[iout:end];contour[1:iout-1]]
    end
end
    
"""
Check `P` is in countour list if so return the index of the contour that contains `P`, otherwise return false."""
function isincountourlist(P::CartesianIndex, contour_list)
    if typeof(P) <: Tuple
        P = CartesianIndex(P)
    end

    for (i, contour) in enumerate(contour_list)
        if P in contour
            return i
        end
    end
    
    return false
    
end

# finds direction between two given pixels
function from_to(from, to, dir_delta)
    delta = to-from
    return findall(x->x == delta, dir_delta)[1]
end

function detect_move(image, p0, p2, nbd, border, done, dir_delta)
    dir = from_to(p0, p2, dir_delta)
    moved = clockwise(dir)
    p1 = CartesianIndex(0, 0)
    while moved != dir ## 3.1
        newp = move(p0, image, moved, dir_delta)
        if newp[1]!=0
            p1 = newp
            break
        end
        moved = clockwise(moved)
    end

    if p1 == CartesianIndex(0, 0)
        return
    end

    p2 = p1 ## 3.2
    p3 = p0 ## 3.2
    done .= false
    while true
        dir = from_to(p3, p2, dir_delta)
        moved = counterclockwise(dir)
        p4 = CartesianIndex(0, 0)
        done .= false
        while true ## 3.3
            p4 = move(p3, image, moved, dir_delta)
            if p4[1] != 0
                break
            end
            done[moved] = true
            moved = counterclockwise(moved)
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

function clockwise(dir)
    return (dir)%8 + 1
end


function counterclockwise(dir)
    return (dir+6)%8 + 1
end

# move from current pixel to next in given direction
function move(pixel, image, dir, dir_delta)
    newp = pixel + dir_delta[dir]
    height, width = size(image)
    if (0 < newp[1] <= height) &&  (0 < newp[2] <= width)
        if image[newp]!=0
            return newp
        end
    end
    return CartesianIndex(0, 0)
end
