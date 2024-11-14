function watershed1(bw::T) where {T<:Union{BitMatrix,AbstractMatrix{Bool}}}
    seg = -IceFloeTracker.bwdist(.!bw)
    mask2 = imextendedmin(seg)
    seg = impose_minima(seg, mask2)
    cc = label_components(imregionalmin(seg), trues(3, 3))
    w = ImageSegmentation.watershed(seg, cc)
    lmap = labels_map(w)
    return Images.isboundary(lmap)
end

function _reconst_watershed(morph_residue::Matrix{<:Integer}, se::Matrix{Bool}=se_disk20())
    mr_reconst = to_uint8(IceFloeTracker.reconstruct(morph_residue, se, "erosion", false))
    mr_reconst .= to_uint8(IceFloeTracker.reconstruct(mr_reconst, se, "dilation", true))
    mr_reconst .= imcomplement(mr_reconst)
    return mr_reconst
end

function watershed2(morph_residue, segment_mask, ice_mask)
    # Task 1: Reconstruct morph_residue
    task1 = Threads.@spawn begin
        mr_reconst = _reconst_watershed(morph_residue)
        mr_reconst = ImageMorphology.local_maxima(mr_reconst; connectivity=2) .> 0
    end

    # Task 2: Calculate gradient magnitude
    task2 = Threads.@spawn begin
        gmag = imgradientmag(histeq(morph_residue))
    end

    # Wait for both tasks to complete
    mr_reconst = fetch(task1)
    gmag = fetch(task2)

    minimamarkers = Bool.(mr_reconst) .| segment_mask .| ice_mask
    gmag .= impose_minima(gmag, minimamarkers)
    cc = label_components(imregionalmin(gmag), trues(3, 3))
    w = ImageSegmentation.watershed(morph_residue, cc)
    lmap = labels_map(w)
    return (fgm=mr_reconst, L0mask=isboundary(lmap))
end
