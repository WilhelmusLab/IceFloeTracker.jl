import Images: erode, local_minima, mreconstruct

"""
    imhmin(img, h)

Mimics MATLAB's imhmin function that computes the H-minima transform, which suppresses all regional minima of `img` whose depth is less than `h`. Computed as the reconstruction by erosion of `img .+ h` above `img`, with 8-connectivity in 2D.

# Arguments
- `img`: grayscale image
- `h`: suppress minima below this depth threshold

# Examples
```jldoctest
julia> img = fill(10.0, 5, 5); img[2:4, 2:4] .= [7 7 7; 7 2 7; 7 7 7];

julia> imhmin(img, 4)
5×5 Matrix{Float64}:
 10.0  10.0  10.0  10.0  10.0
 10.0   7.0   7.0   7.0  10.0
 10.0   7.0   6.0   7.0  10.0
 10.0   7.0   7.0   7.0  10.0
 10.0  10.0  10.0  10.0  10.0
```
"""
function imhmin(img, h)
    raised = img .+ h
    return mreconstruct(erode, raised, img)
end

"""
    imextendedmin(img)

Mimics MATLAB's imextendedmin function that computes the extended-minima transform, which is the regional minima of the H-minima transform. Regional minima are connected components of pixels with a constant intensity value. This function returns a transformed bitmatrix.

# Arguments
- `img`: image object
- `h`: suppress minima below this depth threshold
- `conn`: neighborhood connectivity; in 2D 1 = 4-neighborhood and 2 = 8-neighborhood
"""
function imextendedmin(img::AbstractArray, h::Int=2, conn::Int=2)::BitMatrix
    mask = imhmin(img, h)
    mask_minima = local_minima(mask; connectivity=conn)
    return mask_minima .> 0
end

"""
    imregionalmin(img, conn=2)

Compute the regional minima of the image `img` using the connectivity `conn`.

Returns a bitmatrix of the same size as `img` with the regional minima.

# Arguments
- `img`: Image object
- `conn`: Neighborhood connectivity; in 2D, 1 = 4-neighborhood and 2 = 8-neighborhood
"""
function imregionalmin(img, conn=2)
    return local_minima(img; connectivity=conn) .> 0
end
