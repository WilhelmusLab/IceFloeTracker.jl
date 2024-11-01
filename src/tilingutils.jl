"""
    getfit(dims::Tuple{Int,Int}, side_length::Int)::Tuple{Int,Int}

Calculate how many tiles of a given side length fit into the given dimensions.

# Arguments
- `dims::Tuple{Int,Int}`: A tuple representing the dimensions (width, height).
- `side_length::Int`: The side length of the tile.

# Returns
- `Tuple{Int,Int}`: A tuple representing the number of tiles that fit along each dimension.

# Examples
```
julia> getfit((10, 20), 5)
(2, 4)

julia> getfit((15, 25), 5)
(3, 5)
"""
function getfit(dims::Tuple{Int,Int}, side_length::Int)::Tuple{Int,Int}
    return dims .÷ side_length
end

"""
    get_area_missed(side_length::Int, dims::Tuple{Int,Int})::Float64

Calculate the proportion of the area that is not covered by tiles of a given side length.

# Arguments
- `side_length::Int`: The side length of the tile.
- `dims::Tuple{Int,Int}`: A tuple representing the dimensions (width, height).

# Returns
- `Float64`: The proportion of the area that is not covered by the tiles.

# Examples
```
julia> get_area_missed(5, (10, 20))
0.0

julia> get_area_missed(7, (10, 20))
0.51
"""
function get_area_missed(side_length::Int, dims::Tuple{Int,Int})::Float64
    area = prod(dims)
    return 1 - prod(getfit(dims, side_length)) * side_length^2 / area
end

"""
    get_optimal_tile_size(l0::Int, dims::Tuple{Int,Int}) -> Int

Calculate the optimal tile size in the range [l0-1, l0+1] for the given size `l0` and image dimensions `dims`.

# Description
This function computes the optimal tile size for tiling an area with given dimensions. It ensures that the initial tile size `l0` is at least 2 and not larger than any of the given dimensions. The function evaluates candidate tile sizes and selects the one that minimizes the area missed by its corresponding tiling. In case of a tie, it prefers the larger tile size.

# Example
```
julia> get_optimal_tile_size(3, (10, 7))
2
B```
"""
function get_optimal_tile_size(l0::Int, dims::Tuple{Int,Int})::Int
    l0 < 2 && error("l0 must be at least 2")
    any(l0 .> dims) && error("l0 = $l0 is too large for the given dimensions $dims")

    minimal_shift = l0 == 2 ? 0 : 1
    candidates = [l0 + i for i in (-minimal_shift):1]

    minl, M = 0, Inf
    for side_length in candidates
        missedarea = get_area_missed(side_length, dims)
        if missedarea <= M # prefer larger side_length in case of tie
            M, minl = missedarea, side_length
        end
    end
    return minl
end

"""
    get_tile_meta(tile)

Extracts metadata from a given tile.

# Arguments
- `tile`: A collection of tuples, where each tuple represents a coordinate pair.

# Returns
- A tuple `(a, b, c, d)` where:
  - `a`: The first element of the first tuple in `tile`.
  - `b`: The last element of the first tuple in `tile`.
  - `c`: The first element of the last tuple in `tile`.
  - `d`: The last element of the last tuple in `tile`.
"""
function get_tile_meta(tile)
    a, c = first.(tile)
    b, d = last.(tile)
    return [a, b, c, d]
end

"""
    bump_tile(tile::Tuple{UnitRange{Int64}, UnitRange{Int64}}, dims::Tuple{Int,Int})::Tuple{UnitRange{Int}, UnitRange{Int}}

Adjust the tile dimensions by adding extra rows and columns.

# Arguments
- `tile::Tuple{Int,Int,Int,Int}`: A tuple representing the tile dimensions (a, b, c, d).
- `dims::Tuple{Int,Int}`: A tuple representing the extra rows and columns to add (extrarows, extracols).

# Returns
- `Tuple{UnitRange{Int}, UnitRange{Int}}`: A tuple of ranges representing the new tile dimensions.

# Examples
```julia
julia> bump_tile((1:3, 1:4), (1, 1))
(1:4, 1:5)
"""
function bump_tile(tile::Tuple{UnitRange{S},UnitRange{S}}, dims::Tuple{S,S}) where {S<:Int}
    extrarows, extracols = dims
    a, b, c, d = get_tile_meta(tile)
    b += extrarows
    d += extracols
    return (a:b, c:d)
end

"""
    get_tile_dims(tile)

Calculate the dimensions of a tile.

# Arguments
- `tile::Tuple{UnitRange{Int},UnitRange{Int}}`: A tuple representing the tile dimensions.

# Returns
- `Tuple{Int,Int}`: A tuple representing the width and height of the tile.

# Examples
```julia
julia> get_tile_dims((1:3, 1:4))
(4, 3)
"""
function get_tile_dims(tile)
    a, b, c, d = get_tile_meta(tile)
    width, height = d - c + 1, b - a + 1
    return (width, height)
end

"""
    get_tiles(array, t::Tuple{Int,Int})

Generate a collection of tiles from an array.

The function adjusts the bottom and right edges of the tile matrix if they are smaller than half the tile sizes in `t`.
"""
function get_tiles(array, t::Tuple{T,T}) where {T<:Union{Int,Int64}}
    a, b = t
    tiles = collect(TileIterator(axes(array), (a, b)))
    _a, _b = size(array)

    bottombump = mod(_a, a)
    rightbump = mod(_b, b)

    if bottombump == 0 && rightbump == 0
        return tiles
    end

    crop_height, crop_width = 0, 0

    # Adjust bottom edge if necessary
    if bottombump <= a ÷ 2
        bottom_edge = tiles[end - 1, :]
        tiles[end - 1, :] .= bump_tile.(bottom_edge, Ref((bottombump, 0)))
        crop_height += 1
    end

    # Adjust right edge if necessary
    if rightbump <= b ÷ 2
        right_edge = tiles[:, end - 1]
        tiles[:, end - 1] .= bump_tile.(right_edge, Ref((0, rightbump)))
        crop_width += 1
    end

    return tiles[1:(end - crop_height), 1:(end - crop_width)]
end

"""
    get_tiles(array, side_length)

Generate a collection of tiles from an array.

Unlike `TileIterator`, the function adjusts the bottom and right edges of the tile matrix if they are smaller than half the tile size `side_length`.
"""
function get_tiles(array, side_length::Int)
    return get_tiles(array, (side_length, side_length))
end

"""
    get_tiles(array; rblocks, cblocks)

Generate a collection of tiles from an array.

The function divides the array into `rblocks` rows and `cblocks` columns of tiles.
"""
function get_tiles(array; rblocks, cblocks)
    rtile, ctile = size(array)
    tile_size = (rtile ÷ rblocks, ctile ÷ cblocks)
    return TileIterator(axes(array), tile_size)
end

"""
    get_brighten_mask(equalized_gray_reconstructed_img, gamma_green)

# Arguments

- `equalized_gray_reconstructed_img`: The equalized gray reconstructed image (uint8 in Matlab).
- `gamma_green`: The gamma value for the green channel (also uint8).

# Returns
Difference equalized_gray_reconstructed_img - gamma_green clamped between 0 and 255.

"""
function get_brighten_mask(equalized_gray_reconstructed_img, gamma_green)
    return to_uint8(equalized_gray_reconstructed_img - gamma_green)
end

"""
    imbrighten(img, brighten_mask, bright_factor)

Brighten the image using a mask and a brightening factor.

# Arguments
- `img`: The input image.
- `brighten_mask`: A mask indicating the pixels to brighten.
- `bright_factor`: The factor by which to brighten the pixels.

# Returns
- The brightened image.
"""
function imbrighten(img, brighten_mask, bright_factor)
    img = Float64.(img)
    brighten_mask = brighten_mask .> 0
    img[brighten_mask] .= img[brighten_mask] * bright_factor
    return img = to_uint8(img)
end

function imhist(img, imgtype="uint8")

    # TODO: add validation for arr: either uint8 0:255 or grayscale 0:1

    rng = imgtype == "uint8" ? range(0, 255) : range(0; stop=1, length=256)
    # use range(0, stop=1, length=256) for grayscale images

    # build histogram
    d = Dict(k => 0 for k in rng)
    for i in img
        d[i] = d[i] + 1
    end

    # sort by key (bins)
    k, heights = collect.([Base.keys(d), Base.values(d)])
    order = sortperm(k)
    k, heights = k[order], heights[order]

    return k, heights
end

function get_image_peaks(arr, imgtype="uint8")
    _, heights = imhist(arr, imgtype)

    locs, heights, _ = Peaks.findmaxima(heights)

    # TODO: make this conditional on input args
    order = sortperm(heights; rev=true)
    locs, heights = locs[order], heights[order]

    return (locs=locs, heights=heights)
end

function get_ice_labels_mask(ref_img::Matrix{RGB{N0f8}}, thresholds, factor=1)
    cv = channelview(ref_img)
    cv = [float64.(cv[i, :, :]) .* factor for i in 1:3]
    mask_ice_band_7 = cv[1] .< thresholds[1]
    mask_ice_band_2 = cv[2] .> thresholds[2]
    mask_ice_band_1 = cv[3] .> thresholds[3]
    mask = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
    @debug "Found $(sum(mask)) ice pixels"
    return mask
end

# TODO: add test
function get_nlabel(
    ref_img,
    morph_residue_labels,
    tile,
    factor;
    band_7_threshold::T=5,
    band_2_threshold::T=230,
    band_1_threshold::T=240,
    band_7_threshold_relaxed::T=10,
    band_1_threshold_relaxed::T=190,
    possible_ice_threshold::T=75,
) where {T<:Integer}
    _getnlabel(morphr, tile, mask) = StatsBase.mode(morphr[tile...][mask])

    # Initial attempt to get ice labels
    thresholds = (band_7_threshold, band_2_threshold, band_1_threshold)
    ice_labels_mask = get_ice_labels_mask(ref_img[tile...], thresholds, 255)
    sum(ice_labels_mask) > 1 &&
        return _getnlabel(morph_residue_labels, tile, ice_labels_mask)

    # First relaxation
    thresholds = (band_7_threshold_relaxed, band_2_threshold, band_1_threshold_relaxed)
    ice_labels_mask = get_ice_labels_mask(ref_img[tile...], thresholds, 255)
    sum(ice_labels_mask) > 0 &&
        return _getnlabel(morph_residue_labels, tile, ice_labels_mask)

    # Second/Third relaxation
    return get_nlabel_relaxation(
        ref_img,
        morph_residue_labels,
        tile,
        factor,
        possible_ice_threshold,
        band_7_threshold_relaxed,
        band_2_threshold,
    )
end

function get_nlabel_relaxation(
    ref_img,
    morph_residue_labels,
    tile,
    factor,
    possible_ice_threshold,
    band_7_threshold_relaxed,
    band_2_threshold,
)
    ref_img = ref_img[tile...]
    morph_residue_labels = morph_residue_labels[tile...]

    # filter b/c channels (landmasked channels 2 and 3) and compute peaks
    b, c = [float64.(channelview(ref_img)[i, :, :]) .* factor for i in 2:3]
    b[b .< possible_ice_threshold] .= 0
    c[c .< possible_ice_threshold] .= 0
    pksb, pksc = get_image_peaks.([b, c])

    # return early if no peaks are found
    !all(length.([pksb.locs, pksc.locs]) .> 2) && return 1

    relaxed_thresholds = [band_7_threshold_relaxed, pksb.locs[2], pksc.locs[2]]
    ice_labels = get_ice_labels_mask(ref_img[tile...], relaxed_thresholds, factor)

    sum(ice_labels) > 0 && return StatsBase.mode(morph_residue_labels[ice_labels])

    # Final relaxation
    mask_b = b .> band_2_threshold
    sum(mask_b) > 0 && return StatsBase.mode(morph_residue_labels[mask_b])

    # TODO: Should a fallback value be added? Return nothing if no ice is found? return 1? throw error?
end

function watershed(bw::T) where {T<:Union{BitMatrix,AbstractMatrix{Bool}}}
    seg = -IceFloeTracker.bwdist(.!bw)
    mask2 = imextendedmin(seg, 2)
    seg = impose_minima(seg, mask2)
    cc = label_components(imregionalmin(seg), trues(3, 3))
    w = ImageSegmentation.watershed(seg, cc)
    lmap = labels_map(w)
    return Images.isboundary(lmap)
end

function imextendedmin(img, h)
    return imregionalmin(ImageSegmentation.hmin_transform(img, h))
end

function imregionalmin(A, conn=2)
    return ImageMorphology.local_minima(A; connectivity=conn) .> 0
end

function impose_minima(I::AbstractArray{T}, BW::AbstractArray{Bool}) where {T<:Integer}
    marker = 255 .* BW
    mask = imcomplement(min.(I .+ 1, 255 .- marker))
    reconstructed = IceFloeTracker.MorphSE.mreconstruct(
        IceFloeTracker.MorphSE.dilate, marker, mask
    )
    return IceFloeTracker.imcomplement(Int.(reconstructed))
end

function impose_minima(
    I::AbstractArray{T}, BW::AbstractMatrix{Bool}
) where {T<:AbstractFloat}
    # compute shift
    b, a = extrema(I)
    rng = b - a
    h = rng == 0 ? 0.1 : rng / 1000

    marker = -Inf * BW .+ (Inf * .!BW)
    mask = min.(I .+ h, marker)

    return 1 .- IceFloeTracker.MorphSE.mreconstruct(
        IceFloeTracker.MorphSE.dilate, 1 .- marker, 1 .- mask
    )
end

# TODO: Add tests for get_new2 and get_new3
"""
    getnew2(morph_residue, local_maxima_mask, factor, segment_mask, L0mask)


Calculate the new image `new2` from the input image `morph_residue`.

# Arguments
- `morph_residue`: The morphological residue image.
- `local_maxima_mask`: The local maxima mask.
- `factor`: The factor to apply to the local maxima mask.
- `segment_mask`: The segment mask -- intersection of bw1 and bw2 in first tiled workflow of `master.m`.
- `L0mask`: zero-labeled pixels from watershed.
"""
function get_new2(morph_residue, local_maxima_mask, factor, segment_mask, L0mask)
    new2 = to_uint8(morph_residue .+ local_maxima_mask .* factor)
    new2[segment_mask .|| L0mask] .= 0
    return MorphSE.fill_holes(new2)
end

"""
    get_new3(new2, L0mask, radius, amount, local_maxima_mask, factor, segment_mask)

Calculate the new image `new3` from the input image `new2`.

# Arguments
- `img`: The input image.
- `L0mask`: zero-labeled pixels from watershed.
- `radius`: The radius of the unsharp mask.
- `amount`: The amount of unsharp mask.
- `local_maxima_mask`: The local maxima mask.
- `factor`: The factor to apply to the local maxima mask.
- `segment_mask`: The segment mask -- intersection of bw1 and bw2 in first tiled workflow of `master.m`.

"""
function get_new3(img, L0mask, radius, amount, local_maxima_mask, factor, segment_mask)
    new3 = unsharp_mask(img, radius, amount, 255)
    new3[L0mask] .= 0
    new3 = reconstruct(new3, se, "dilation", false)
    new3[segment_mask] .= 0
    return to_uint8(new3 + local_maxima_mask .* factor)
end

"""
    get_final(img, label, segment_mask, se_erosion, se_dilation)

Final processing following the tiling workflow.
"""
function get_final(img, label, segment_mask, se_erosion, se_dilation)
    img = hbreak(img)

    # doesn't do much/anything
    if label == 1
        img[1:1] .= false
    end

    # slow for big images
    img = morph_fill(img)

    # only works for label 1, whose value tends to be arbirary
    if label == 1
        img[segment_mask] .= false
    end

    # tends to fill more than matlabs imfill
    img = IceFloeTracker.MorphSE.fill_holes(img)

    marker = branch(img)

    mask = MorphSE.erode(marker, se_erosion)
    mask = MorphSE.dilate(mask, se_dilation)

    # doesnt do much
    if label == 1
        mask[1] = false
    end

    final = MorphSE.mreconstruct(MorphSE.dilate, marker, mask)
    return final
end
