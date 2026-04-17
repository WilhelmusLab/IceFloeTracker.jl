import Images: centered
"""
    _generate_se!(se)

Generate a structuring element by leveraging symmetry (mirroring and inverting) a given initial structuring element.
"""
function _generate_se!(se)
    for d in [1, 2]
        se .= se .| reverse(se; dims=d)
    end
    se .= .!se
    return nothing
end

"""
    se_octagon(r)

Construct an octagonal structuring element with radius r.

""" # dmw: should this have the same offset indexing we get in ImageMorphology?
function strel_octagon(r)
    se = [sum(c.I) <= r/2 + 2 for c in CartesianIndices((2*r + 1, 2*r + 1))]
    _generate_se!(se)
    return centered(se)
end
"""
    strel_disk(r)

Construct a disk-shaped structuring element with radius r, diameter 2r+1.
"""
function strel_disk(r)
    se = [sum(abs.(c.I .- (r + 1)) .^ 2) for c in CartesianIndices((2*r + 1, 2*r + 1))]
    return centered(sqrt.(se) .<= r)
end
