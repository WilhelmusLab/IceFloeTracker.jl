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

export imregionalmin

end
