function _generate_se(se)
    se = se .| reverse(se; dims=1)
    se = se .| reverse(se; dims=2)
    return .!se
end

make_landmask_se() = begin
    se = [sum(c.I) <= 29 for c in CartesianIndices((99, 99))]
    _generate_se(se)
end

se_disk20 = make_landmask_se

se_disk4() = begin
    se = zeros(Bool, 7, 7)
    se[4, 4] = 1
    bwdist(se) .<= 3.6
end

se_disk20() = begin
    se = [sum(c.I) < 11 for c in CartesianIndices((39, 39))]
    return _generate_se(se)
end
