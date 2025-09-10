function watershed1(bw::T) where {T<:Union{BitMatrix,AbstractMatrix{Bool}}}
    seg = -IceFloeTracker.bwdist(.!bw)
    mask2 = imextendedmin(seg)
    seg = impose_minima(seg, mask2)
    cc = label_components(imregionalmin(seg), trues(3, 3))
    w = ImageSegmentation.watershed(seg, cc)
    lmap = labels_map(w)
    return Images.isboundary(lmap) .> 0
end


# dmw: Expose this function in the watershed workflow directly, using ImageMorphology functions
function _reconst_watershed(morph_residue::Matrix{<:Integer}, se::Matrix{Bool}=se_disk20())
    mr_reconst = to_uint8(IceFloeTracker.reconstruct(morph_residue, se, "erosion", false))
    mr_reconst .= to_uint8(IceFloeTracker.reconstruct(mr_reconst, se, "dilation", true))
    mr_reconst .= imcomplement(mr_reconst)
    return mr_reconst
end

# dmw: figure out new name for this that is clearer that "watershed2".
function watershed2(morph_residue, segment_mask, ice_mask)
    # TODO: reconfigure to use async tasks or threads
    # Task 1: Reconstruct morph_residue
    # task1 = Threads.@spawn begin

    # dmw: temporary fix; TBD rewrite _reconst_watershed to use ImageMorphology directly
    maximum(morph_residue) <= 1. && (
        morph_residue = Int64.(round.(Float64.(morph_residue) .* 255, digits=0))
    )

    # dmw: converting this to binary makes the meaning different, it's now a binary mask
    mr_reconst = _reconst_watershed(morph_residue)
    mr_reconst = ImageMorphology.local_maxima(mr_reconst; connectivity=2) .> 0
    # end

    # Task 2: Calculate gradient magnitude
    # task2 = Threads.@spawn begin
    gmag = imgradientmag(histeq(morph_residue))
    # end

    # Wait for both tasks to complete
    # mr_reconst = fetch(task1)
    # gmag = fetch(task2)

    minimamarkers = mr_reconst .| segment_mask .| ice_mask
    gmag .= impose_minima(gmag, minimamarkers)
    cc = label_components(imregionalmin(gmag), trues(3, 3))
    w = ImageSegmentation.watershed(morph_residue, cc)
    lmap = labels_map(w) # dmw: inspect results here. Should we not multiply by an earlier mask to remove weird edges?
    return (fgm=mr_reconst, L0mask=isboundary(lmap) .> 0) # dmw: what does fgm stand for?
end
