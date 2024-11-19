"""
    _generate_se!(se)

Generate a structuring element by leveraging symmetry (mirroring and inverting) a given initial structuring element.
"""
function _generate_se!(se)
    se .= se .| reverse(se; dims=1)
    se .= se .| reverse(se; dims=2)
    se .= .!se
    return nothing
end

function se_disk50()
    se = [sum(c.I) <= 29 for c in CartesianIndices((99, 99))]
    _generate_se!(se)
    return se
end

make_landmask_se = se_disk50

function se_disk4()
    se = zeros(Bool, 7, 7)
    se[4, 4] = 1
    return bwdist(se) .<= 3.6
end

function se_disk20()
    se = [sum(c.I) <= 11 for c in CartesianIndices((39, 39))]
    _generate_se!(se)
    return se
end

function se_disk2()
    se = [sum(c.I) <= 3 for c in CartesianIndices((5,5))]
    _generate_se!(se)
    return se
end
