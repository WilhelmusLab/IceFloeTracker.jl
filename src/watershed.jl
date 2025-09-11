function watershed1(bw::T) where {T<:Union{BitMatrix,AbstractMatrix{Bool}}}
    seg = -IceFloeTracker.bwdist(.!bw)
    mask2 = imextendedmin(seg)
    seg = impose_minima(seg, mask2)
    cc = label_components(imregionalmin(seg), trues(3, 3))
    w = ImageSegmentation.watershed(seg, cc)
    lmap = labels_map(w)
    return Images.isboundary(lmap) .> 0
    # dmw: We actually could get a really nice segmented image from returning the segmentation
    # defined by labels_map(w) .* bw; this keeps the segments separated by watershed without imposing
    # minimum separation.
end




# dmw: figure out new name for this that is clearer that "watershed2".
function watershed2(morph_residue, segment_mask, ice_mask)

    # Two passes of reconstruction by dilation, first eroding to find markers, then dilating the complement.
    # This is a really large structuring element, and the erosion of the morphological residue
    # is very similar to the dilation of the complement.
    mr_reconst = mreconstruct(dilate,
                              erode(morph_residue, se_disk20()),
                              morph_residue,
                              collect(strel_box((3, 3)))
                            ) # 3x3 box is default for sk morphology
    
    # dmw: question with local max: are we trying to find places between floes, or find the floes themselves?
    # Originally the function was returning the complement of the reconstruction of the complement:
    # i.e., inverting the image, reconstructing it, then inverting it back.
    # However, then the maxima are floes, not gaps.
    # So here I've taken out the final complement step.
    mr_reconst .= complement.(
                        mreconstruct(dilate,
                                     complement.(dilate(mr_reconst, se_disk20())), 
                                     complement.(mr_reconst),
                                     collect(strel_box((3, 3)))
                                    )
                            )


    # dmw: The local_maxima function selects the level sets larger than their surroundings.
    # It should be identifying ice floes, but may leave many ice floes out.
    local_max = ImageMorphology.local_maxima(mr_reconst; connectivity=2) .> 0
    
    
    # Compute the image gradients of the equalized morphological residue using the Sobel operator
    mr_equalized = adjust_histogram(morph_residue, Equalization(nbins = 256, minval = 0, maxval = 1))
    Gx, Gy = imgradients(mr_equalized,  KernelFactors.sobel, "replicate")
    gmag = hypot.(Gx, Gy)
    gmag .= gmag ./ maximum(gmag)


    # dmw: this part doesn't make sense to me. Why would we want the borders (positive) and the local maxima and the ice mask?
    # segment_mask .|| ice_mask would be almost all 1s.
    # What is expected from the ice mask to make this step make sense?
    minimamarkers = local_max .| segment_mask .| ice_mask
    gmag .= impose_minima(gmag, minimamarkers)
    cc = label_components(imregionalmin(gmag), trues(3, 3))
    w = ImageSegmentation.watershed(morph_residue, cc)
    lmap = labels_map(w) # dmw: inspect results here. Should we not multiply by an earlier mask to remove weird edges?
    return (fgm=local_max, L0mask=isboundary(lmap) .> 0) # dmw: what does fgm stand for?
end
