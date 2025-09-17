# Helper functions

"""
    loadimg(; dir::String, fname::String)

Load an image from `dir` with filename `fname` into a matrix of `Float64` values. Returns the loaded image.
"""
function loadimg(; dir::String, fname::String)
    return (x -> float64.(x))(load(joinpath(dir, fname)))
end

"""
    add_padding(img, style)

Extrapolate the image `img` according to the `style` specifications type. Returns the extrapolated image.

# Arguments
- `img`: Image to be padded.
- `style`: A supported type (such as `Pad` or `Fill`) representing the extrapolation style. See the relevant [documentation](https://juliaimages.org/latest/function_reference/#ImageFiltering) for details.

See also [`remove_padding`](@ref)
"""
function add_padding(img, style::Union{Pad,Fill})::Matrix
    return collect(Images.padarray(img, style))
end

"""
    remove_padding(paddedimg, border_spec)

Removes padding from the boundary of padded image `paddedimg` according to the border specification `border_spec` type. Returns the cropped image.

# Arguments
- `paddedimg`: Pre-padded image.
- `border_spec`: Type representing the style of padding (such as `Pad` or `Fill`) with which `paddedimg` is assumend to be pre-padded. Example: `Pad((1,2), (3,4))` specifies 1 row on the top, 2 columns on the left, 3 rows on the bottom, and 4 columns on the right boundary.

See also [`add_padding`](@ref)
"""
function remove_padding(paddedimg, border_spec::Union{Pad,Fill})::Matrix
    top, left = border_spec.lo
    bottom, right = border_spec.hi
    return paddedimg[(top + 1):(end - bottom), (left + 1):(end - right)]
end

"""
    impose_minima(I::AbstractArray{T}, BW::AbstractArray{Bool}) where {T<:Integer}

Use morphological reconstruction to enforce minima on the input image `I` at the positions where the binary mask `BW` is non-zero.

It supports both integer and grayscale images using different implementations for each.
"""
function impose_minima(I::AbstractArray{T}, BW::AbstractArray{Bool}) where {T<:Integer}
    marker = 255 .* BW
    mask = imcomplement(min.(I .+ 1, 255 .- marker))
    reconstructed = sk_morphology.reconstruction(marker, mask)
    return IceFloeTracker.imcomplement(Int.(reconstructed))
end

function impose_minima(
    I::AbstractArray{T}, BW::AbstractMatrix{Bool}
) where {T<:AbstractFloat}
    # compute shift
    a, b = extrema(I)
    rng = b - a
    h = rng == 0 ? 0.1 : rng / 1000

    marker = -Inf * BW .+ (Inf * .!BW)
    mask = min.(I .+ h, marker)

    return 1 .- sk_morphology.reconstruction(1 .- marker, 1 .- mask)
end

"""
    padnhood(img, I, nhood)

Pad the matrix `img[nhood]` with zeros according to the position of `I` within the edges`img`.

Returns `img[nhood]` if `I` is not an edge index.
"""
function padnhood(img, I, nhood)
    # adaptive padding
    maxr, maxc = size(img)
    tofill = SizedMatrix{3,3}(zeros(Int, 3, 3))
    @views if I == CartesianIndex(1, 1) # top left corner`
        tofill[2:3, 2:3] = img[nhood]
    elseif I == CartesianIndex(maxr, 1) # bottom left corner
        tofill[1:2, 2:3] = img[nhood]
    elseif I == CartesianIndex(1, maxc) # top right corner
        tofill[2:3, 1:2] = img[nhood]
    elseif I == CartesianIndex(maxr, maxc) # bottom right corner
        tofill[1:2, 1:2] = img[nhood]
    elseif I[1] == 1 # top edge (first row)
        tofill[2:3, 1:3] = img[nhood]
    elseif I[2] == 1 # left edge (first col)
        tofill[1:3, 2:3] = img[nhood]
    elseif I[1] == maxr # bottom edge (last row)
        tofill[1:2, 1:3] = img[nhood]
    elseif I[2] == maxc # right edge (last row)
        tofill[1:3, 1:2] = img[nhood]
    else
        tofill = img[nhood]
    end
    return tofill
end

"""
    _bin9todec(v)

Get decimal representation of a bit vector `v` with the leading bit at its leftmost posistion.

Example
```
julia> _bin9todec([0 0 0 0 0 0 0 0 0])
0

julia> _bin9todec([1 1 1 1 1 1 1 1 1])
511
```
"""
function _bin9todec(v::AbstractArray)::Int64
    return sum(vec(v) .* 2 .^ (0:(length(v) - 1)))
end

"""
    _operator_lut(I, img, nhood, lut1, lut2)

Look up the neighborhood `nhood` in lookup tables `lut1` and `lut2`.

Handles cases when the center of `nhood` is on the edge of `img` using data in `I`.
"""
function _operator_lut(
    I::CartesianIndex{2},
    img::AbstractArray{Bool},
    nhood::CartesianIndices{2,Tuple{UnitRange{Int64},UnitRange{Int64}}},
    lut1::Vector{Int64},
    lut2::Vector{Int64},
)::SVector{2,Int64}

    # corner pixels
    length(nhood) == 4 && return @SVector [false, 0]

    val = IceFloeTracker._bin9todec(_pad_handler(I, img, nhood)) + 1

    return @SVector [lut1[val], lut2[val]]
end

function _operator_lut(
    I::CartesianIndex{2},
    img::AbstractArray{Bool},
    nhood::CartesianIndices{2,Tuple{UnitRange{Int64},UnitRange{Int64}}},
    lut::Vector{T},
)::T where {T} # for bridge

    # corner pixels
    length(nhood) == 4 && return false # for bridge and some other operations like hbreak, branch

    return lut[_bin9todec(_pad_handler(I, img, nhood)) + 1]
end

function _pad_handler(I, img, nhood)
    (length(nhood) == 6) && return padnhood(img, I, nhood) # edge pixels
    return @view img[nhood]
end

"""
    callable_store()

Create a store and a callback function to add key-value pairs to the store.

Returns a `store::Dict` and a `callback::Function` which stores any kwargs passed to it in the `store`.

# Examples

Basic usage is to store values using the callback function
```julia-repl
julia> store, callback = callable_store()
julia> store
Dict{Any, Any}()

julia> callback(;foo="bar")  # echoes the updated store
Dict{Any, Any} with 1 entry:
  :foo => "bar"

julia> store  # values are available from the store object
Dict{Any, Any} with 1 entry:
  :foo => "bar"
```

A real-world use case is to extract data from a segmentation algorithm run:
```julia-repl
julia> intermediate_results, intermediate_results_callback = callable_store()
julia> data = first(Watkins2025GitHub(; ref="a451cd5e62a10309a9640fbbe6b32a236fcebc70")());
julia> segments = LopezAcosta2019Tiling()(
    data.modis_truecolor,
    data.modis_falsecolor,
    data.modis_landmask;
    intermediate_results_callback,
)
Segmented Image with:
  labels map: 400×400 Matrix{Int64}
  number of labels: 12

julia> intermediate_data
Dict{Any, Any} with 16 entries:
  :binarized_tiling                       => Bool[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :icemask                                => Bool[1 1 … 1 1; 1 1 … 1 1; … ; 0 0 … 1 1; 0 0 … 1 1]
  :equalized_gray                         => [0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :morphed_residue                        => [0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :L0mask                                 => Bool[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :segmented                              => Segmented Image with:…
  :prelim_icemask2                        => [255 255 … 255 255; 255 255 … 255 255; … ; 255 255 … 255 255; 255 255 … 255 255]
  :equalized_gray_sharpened_reconstructed => [0 0 … 0 0; 0 0 … 0 0; … ; 255 255 … 255 255; 255 255 … 255 255]
  :gammagreen                             => [190.35 190.23 … 182.93 185.03; 191.68 190.6 … 185.04 192.08; … ; 163.87 173.33 … 108.02 108.18; 166.14 173.3 … 112.35 112.32]
  :segment_mask                           => Bool[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :ref_img_cloudmasked                    => RGB{N0f8}[RGB{N0f8}(0.0,0.0,0.0) RGB{N0f8}(0.0,0.0,0.0) … RGB{N0f8}(0.008,0.706,0.761) RGB{N0f8}(0.0,0.722,0.769); RGB{N0f8}(0.0,0.0,0.0) RGB{N0f8}(0.0,0.0,0.0) … RGB{N0f8}(0.039,0.702,0.784) RGB{N0f8}(0.075,0.784,0.859); … ; RGB{…
  :prelim_icemask                         => Bool[0 0 … 0 0; 0 0 … 0 0; … ; 0 0 … 0 0; 0 0 … 0 0]
  :equalized_gray_reconstructed           => [0 0 … 0 0; 0 0 … 0 0; … ; 255 255 … 255 255; 255 255 … 255 255]
  :final                                  => Bool[0 0 … 0 0; 0 1 … 1 0; … ; 0 0 … 1 0; 0 0 … 0 0]
  :local_maxima_mask                      => [255 255 … 255 255; 255 255 … 255 255; … ; 255 255 … 255 255; 255 255 … 255 255]
  :labeled                                => [0 0 … 0 0; 0 1 … 1 0; … ; 0 0 … 9 0; 0 0 … 0 0]
```
"""
function callable_store()::Tuple{Dict,Function}
    store = Dict()
    function callback(; kwargs...)
        return merge!(store, Dict(kwargs))
    end
    return store, callback
end
