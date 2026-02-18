module CLAHE

using Images: YIQ, channelview, padarray, Pad
using Images.ImageContrastAdjustment:
    AbstractHistogramAdjustmentAlgorithm, GenericGrayImage, imresize, build_histogram, adjust_histogram, adjust_histogram!
using Images.ImageCore
using TiledIteration: TileIterator

"""
```
    ContrastLimitedAdaptiveHistogramEqualization <: AbstractHistogramAdjustmentAlgorithm
    ContrastLimitedAdaptiveHistogramEqualization(; nbins = 128, minval = 0, maxval = 1, rblocks = 8, cblocks = 8, clip = 0.1)

    adjust_histogram([T,] img, f::ContrastLimitedAdaptiveHistogramEqualization)
    adjust_histogram!([out,] img, f::ContrastLimitedAdaptiveHistogramEqualization)
```

Performs Contrast Limited Adaptive Histogram Equalisation (CLAHE) on the input
image.

This version is based on the code in:
GraphicsGems IV, "Contrast Limited Adaptive Histogram Equalization".

# References
1. P. Heckbert (ed.). *Graphics gems IV*. Cambridge MA, USA: AP Professional, 1994.

"""
@kwdef struct ContrastLimitedAdaptiveHistogramEqualization{
    T₁<:Union{Real,AbstractGray},T₂<:Union{Real,AbstractGray},T₃<:Real
} <: AbstractHistogramAdjustmentAlgorithm
    nbins::Int = 128
    minval::T₁ = 0.0
    maxval::T₂ = 1.0
    rblocks::Int = 8
    cblocks::Int = 8
    clip::T₃ = 0.1
end

function (f::ContrastLimitedAdaptiveHistogramEqualization)(
    out::GenericGrayImage, img::GenericGrayImage
)
    validate_parameters(f)
    height, width = length.(axes(img))

    # If necessary, resize the image so that the requested number of blocks fit exactly.
    resized_height = ceil(Int, height / (2 * f.rblocks)) * 2 * f.rblocks
    resized_width = ceil(Int, width / (2 * f.cblocks)) * 2 * f.cblocks
    must_pad = (resized_height != height) || (resized_width != width)
    if must_pad
        left, right = ceil(Int, (resized_width - width) / 2), floor(Int, (resized_width - width) / 2)
        top, bottom = ceil(Int, (resized_height - height) / 2), floor(Int, (resized_height - height) / 2)
        img_padded = padarray(img, Pad(:reflect, (top, left), (bottom, right)))
    else
        img_padded = img
    end
    out_padded = similar(img_padded)

    # Size of each contextual region
    rsize = resized_height ÷ f.rblocks
    csize = resized_width ÷ f.cblocks

    # Calculate actual clip limit
    clip_limit = f.clip * (rsize * csize) / f.nbins
    clip_limit < 1 && (clip_limit = 1)
    clip_limit = Int(floor(clip_limit))

    # Process each contextual region
    histograms = Array{Any}(undef, f.rblocks, f.cblocks)
    tiles = collect(TileIterator(axes(img_padded), (rsize, csize)))
    for (I, tile) in zip(CartesianIndices(tiles), tiles)
        rblock, cblock = Tuple(I)
        region = img_padded[tile...]
        edges, raw_counts = build_histogram(
            region, f.nbins; minval=f.minval, maxval=f.maxval
        )
        redistributed_counts = redistribute_histogram(raw_counts, clip_limit)
        mapping_function = map_histogram(
            edges::AbstractArray, redistributed_counts, f.minval, f.maxval
        )
        histograms[rblock, cblock] = mapping_function
    end

    # Interpolate pixel values
    for rblock in 1:(f.rblocks + 1), cblock in 1:(f.cblocks + 1)
        
        # Get the histograms for each block
        if rblock == 1
            idUr, idBr = 1, 1
        elseif rblock == f.rblocks + 1
            idUr, idBr = f.rblocks, f.rblocks
        else
            idUr, idBr = rblock - 1, rblock
        end
        if cblock == 1
            idLc, idRc = 1, 1
        elseif cblock == f.cblocks + 1
            idLc, idRc = f.cblocks, f.cblocks
        else
            idLc, idRc = cblock - 1, cblock
        end
        histUL, histUR = histograms[idUr, idLc], histograms[idUr, idRc]
        histBL, histBR = histograms[idBr, idLc], histograms[idBr, idRc]

        # Get the size of the block, which may be half the normal size if we're on the edge of the image
        rblockpix = rblock ∈ [1, f.rblocks + 1] ? rsize / 2 : rsize
        rblockoffset = rblock == 1 ? 0 : -rsize / 2
        cblockpix = cblock ∈ [1, f.cblocks + 1] ? csize / 2 : csize
        cblockoffset = cblock == 1 ? 0 : -csize / 2
        
        # Get the block itself, indexing from the origin of the axes of the image (which may be negative if we had to pad)
        rorigin, corigin = minimum.(axes(img_padded))
        rstart = Int((rblock - 1) * rsize + rblockoffset) + rorigin
        rend = Int(rstart + rblockpix) - 1
        cstart = Int((cblock - 1) * csize + cblockoffset) + corigin
        cend = Int(cstart + cblockpix) - 1
        
        region = view(img_padded, rstart:rend, cstart:cend)
        out_region = view(out_padded, rstart:rend, cstart:cend)

        # Calculate new values using each of the histograms from the four surrounding blocks
        resultUL = histUL.(region)
        resultUR = histUR.(region)
        resultBL = histBL.(region)
        resultBR = histBR.(region)

        # Interpolate the results bilinearly
        r₁, r₂ = rstart, rend
        c₁, c₂ = cstart, cend
        r = Array(range(r₁, r₂))
        c = Array(range(c₁, c₂))'

        w₁₁ = @. ((r₂ - r) * (c₂ - c))
        w₁₂ = @. ((r₂ - r) * (c - c₁))
        w₂₁ = @. ((r - r₁) * (c₂ - c))
        w₂₂ = @. ((r - r₁) * (c - c₁))
        wₙ = ((r₂ - r₁) * (c₂ - c₁))

        # Write to the correct region in the output image
        @. out_region =
            (w₁₁ * resultUL + w₁₂ * resultUR + w₂₁ * resultBL + w₂₂ * resultBR) / wₙ
    end

    out .= must_pad ? out_padded[1:height, 1:width] : out_padded
    return out
end

function validate_parameters(f::ContrastLimitedAdaptiveHistogramEqualization)
    !(1 <= f.rblocks && 1 <= f.cblocks) &&
        throw(ArgumentError("At least 1 contextual regions required (1x1 or greater)."))
    return nothing
end

function redistribute_histogram(
    counts::AbstractVector{T}, clip_limit::Int
) where {T<:Integer}
    n_excess = sum(max(0, count - clip_limit) for count in counts)
    n_bins = length(counts)
    n_excess == 0 && return counts
    increment, remainder = divrem(n_excess, n_bins)
    new_counts = @. min(counts, clip_limit) + increment
    for i in 1:remainder
        new_counts[i] += 1
    end
    return new_counts
end

function map_histogram(
    edges::AbstractArray,
    counts::AbstractVector{T},
    minval::Union{Real,AbstractGray},
    maxval::Union{Real,AbstractGray},
) where {T<:Integer}
    n_pixels = sum(counts)
    scale = (maxval - minval) / n_pixels
    mapping = similar(counts, Float64) # Can I make this the same type as minval/maxval?
    cumulative = 0
    for i in eachindex(counts)
        cumulative += counts[i]
        mapping[i] = minval + cumulative * scale
        if mapping[i] > maxval
            mapping[i] = maxval
        end
    end

    function _mapping_function(value)
        index = searchsortedfirst(edges, value)
        return mapping[index - 1] # -1 because mapped_values is an OffsetArray
    end
    return _mapping_function
end

function (f::ContrastLimitedAdaptiveHistogramEqualization)(
    out::AbstractArray{<:Color3}, img::AbstractArray{<:Color3}
)
    T = eltype(img)
    yiq = convert.(YIQ, img)
    yiq_view = channelview(yiq)
    #=
       TODO: Understand the cause and solution of this error.
       When I pass a view I run into this error on Julia 1.1.
       ERROR: ArgumentError: an array of type `Base.ReinterpretArray` shares memory with another argument and must
       make a preventative copy of itself in order to maintain consistent semantics,
       but `copy(A)` returns a new array of type `Array{Float64,3}`. To fix, implement:
       `Base.unaliascopy(A::Base.ReinterpretArray)::typeof(A)`
    =#
    #adjust_histogram!(view(yiq_view,1,:,:), f)
    y = comp1.(yiq)
    adjust_histogram!(y, f)
    yiq_view[1, :, :] .= y
    out .= convert.(T, yiq)
    return out
end

export ContrastLimitedAdaptiveHistogramEqualization, adjust_histogram, adjust_histogram!

end
