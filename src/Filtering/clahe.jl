module CLAHE

using Images
using Images.ImageContrastAdjustment:
    AbstractHistogramAdjustmentAlgorithm, GenericGrayImage, imresize
using Images.ImageCore

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
    must_resize = (resized_height != height) || (resized_width != width)
    if must_resize
        img_tmp = imresize(img, (resized_height, resized_width))
        out_tmp = copy(img_tmp)
    else
        img_tmp = img
        out_tmp = copy(img)
    end

    # Size of each contextual region
    rsize = resized_height / f.rblocks
    csize = resized_width / f.cblocks

    # Calculate actual clip limit
    clip_limit = f.clip * (rsize * csize) / f.nbins
    clip_limit < 1 && (clip_limit = 1)
    clip_limit = Int(floor(clip_limit))

    # Process each contextual region
    histograms = Array{Any}(undef, f.rblocks, f.cblocks)
    for rblock in 1:(f.rblocks), cblock in 1:(f.cblocks)
        rstart = Int((rblock - 1) * rsize) + 1
        rend = Int(rblock * rsize)
        cstart = Int((cblock - 1) * csize) + 1
        cend = Int(cblock * csize)
        region = view(img_tmp, rstart:rend, cstart:cend)
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
    for rblock in 1:(f.rblocks + 1), cblock in 1:(f.cblocks + 1) # use zero-indexing here because we're not in the original format
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

        if rblock ∈ [1, f.rblocks + 1]
            rblockpix = rsize / 2
        else
            rblockpix = rsize
        end

        if rblock == 1
            rblockoffset = 0
        else
            rblockoffset = -(rsize / 2)
        end

        if cblock ∈ [1, f.cblocks + 1]
            cblockpix = csize / 2
        else
            cblockpix = csize
        end

        if cblock == 1
            cblockoffset = 0
        else
            cblockoffset = -(csize / 2)
        end

        histUL, histUR = histograms[idUr, idLc], histograms[idUr, idRc]
        histBL, histBR = histograms[idBr, idLc], histograms[idBr, idRc]

        rstart = Int((rblock - 1) * rsize + rblockoffset) + 1
        rend = Int(rstart + rblockpix) - 1

        cstart = Int((cblock - 1) * csize + cblockoffset) + 1
        cend = Int(cstart + cblockpix) - 1

        region = view(img_tmp, rstart:rend, cstart:cend)
        out_region = view(out_tmp, rstart:rend, cstart:cend)

        resultUL, resultUR = histUL.(region), histUR.(region)
        resultBL, resultBR = histBL.(region), histBR.(region)

        resultUL = histUL.(region)
        resultUR = histUR.(region)
        resultBL = histBL.(region)
        resultBR = histBR.(region)

        x₁, x₂ = rstart, rend
        y₁, y₂ = cstart, cend
        x = Array(range(rstart, rend))
        y = Array(range(cstart, cend))'

        w₁₁ = ((x₂ .- x) .* (y₂ .- y))
        w₁₂ = ((x₂ .- x) .* (y .- y₁))
        w₂₁ = ((x .- x₁) .* (y₂ .- y))
        w₂₂ = ((x .- x₁) .* (y .- y₁))
        wₙ = ((x₂ - x₁) * (y₂ - y₁))

        @. out_region =
            (w₁₁ * resultUL + w₁₂ * resultUR + w₂₁ * resultBL + w₂₂ * resultBR) / wₙ
    end

    out .= must_resize ? imresize(out_tmp, (height, width)) : out_tmp
    return out
end

function validate_parameters(f::ContrastLimitedAdaptiveHistogramEqualization)
    # !(0 <= f.clip <= 1) && throw(ArgumentError("The parameter `clip` must be in the range [0..1]."))
    return !(1 <= f.rblocks && 1 <= f.cblocks) &&
           throw(ArgumentError("At least 1 contextual regions required (1x1 or greater)."))
end

function redistribute_histogram(
    counts::AbstractVector{T}, clip_limit::Int
) where {T<:Integer}
    n_excess = sum(max(0, count - clip_limit) for count in counts)
    n_bins = length(counts)
    if n_excess == 0
        return counts
    end
    increment, remainder = divrem(n_excess, n_bins)
    new_counts = similar(counts)
    for i in eachindex(counts)
        new_counts[i] = min(counts[i], clip_limit) + increment
    end
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
    return out .= convert.(T, yiq)
end

export ContrastLimitedAdaptiveHistogramEqualization, adjust_histogram, adjust_histogram!

end # module CLAHE
