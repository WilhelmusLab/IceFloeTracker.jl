
"""
segmentation_C(segmented_b_filled, segmented_b_ice;)

Joins an intermediate file with the final output of segmentation_B to further contrast potential ice floes, returning a mask of potential ice.

# Arguments
- `segmented_b_filled`: binary cloudmasked and landmasked file `segmentation_b_filled` from `segmentation_b.jl`
- `segmented_b_ice`: binary cloudmasked and landmasked file `segmentation_b_ice` from `segmentation_b.jl`
"""

function segmentation_C(
    segmented_b_filled::BitMatrix, segmented_b_ice::BitMatrix
)::BitMatrix
    segmented_c = segmented_b_ice .* segmented_b_filled

    return segmented_c
end
