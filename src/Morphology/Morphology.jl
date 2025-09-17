module Morphology

import ImageMorphology: ImageMorphology

"""
    imregionalmin(img, conn=2)

Compute the regional minima of the image `img` using the connectivity `conn`.

Returns a bitmatrix of the same size as `img` with the regional minima.

# Arguments
- `img`: Image object
- `conn`: Neighborhood connectivity; in 2D, 1 = 4-neighborhood and 2 = 8-neighborhood
"""
function imregionalmin(img, conn=2)
    return ImageMorphology.local_minima(img; connectivity=conn) .> 0
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
    mask = ImageSegmentation.hmin_transform(img, h)
    mask_minima = Images.local_minima(mask; connectivity=conn)
    return mask_minima .> 0
end

export imregionalmin, imextendedmin

end
