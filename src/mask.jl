using Images, FileIO

"""
    masker(mask::AbstractArray, img::AbstractArray{<:Colorant})
    masker(mask::AbstractArray)
    
Returns a version of `img` with masked pixels made transparent.
If `img` has an alpha channel, it is combined with the mask.

masker(mask) returns a function which can be used to apply the mask.

# Examples

With a BitMatrix type of mask, truthy values in the mask are transparent in the output.
```julia-repl
julia> using Images
julia> hide = true
julia> pass = false
julia> bit_mask = [hide hide pass; pass pass hide]
julia> img = parse.(Colorant, ["red" "green" "blue"; "cyan" "magenta" "yellow"])
julia> masker(bit_mask, img)
```

Using the `masker` in a pipeline is also possible:
```julia-repl
julia> img |> masker(bit_mask)
```

Where the mask is itself an image with transparency, 
areas which are opaque in the mask
are transparent in the output.
This corresopnds to overlaying the mask over the image,
and hiding in the output those areas which were masked.
```julia-repl
julia> hide = AGray(0.5, 1)
julia> pass = AGray(1, 0)
julia> agray_mask = [hide hide pass; pass pass hide]
julia> masker(agray_mask, img)
```

Any pixels which are partially opaque in the mask,
will be partially obscured in the output:
```julia-repl
julia> part = AGray(0.75, 0.5)
julia> partial_mask = [hide part pass; pass part hide]
julia> masker(partial_mask, img)
```

If the image already has transparency, 
this is combined with the mask.
```julia-repl
julia> imga = RGBA.(parse.(Colorant, ["red" "transparent" "blue"; "cyan" "transparent" "yellow"]))
julia> masker(agray_mask, imga)
```

Where the mask is an image without transparency, 
any non-zero pixels are masked:
```julia-repl
julia> gray_mask = Gray.([0.5 0.1 0.0; 0.0 0.0 1.0])
julia> masker(gray_mask, img)
julia> rgb_mask = parse.(Colorant, ["red" "red" "black"; "black" "purple" "green"])
julia> masker(rgb_mask, img)
```

Where the mask is an array of Real numbers,
0-pixels will be completely unmasked,
1-pixels completely masked,
and values between partially masked:
```julia-repl
julia> real_mask = [1.0 0.5 0.1; 0.1 0.2 1.0]
julia> masker(real_mask, img)
```

Where values are outside of the range [0, 1], 
they are clamped to whichever of 0 and 1 is nearer:
```julia-repl
julia> out_of_range_mask = [5 2 0.75; -1 -2 1]
julia> masker(out_of_range_mask, img)
```


"""
function masker(mask::AbstractArray)
    masking_alpha_channel = _mask_to_alpha(mask)
    function _apply(img::AbstractArray{<:Color})
        return alphacolor.(img, masking_alpha_channel)
    end
    function _apply(img::AbstractArray{<:TransparentColor})
        combined_masking_alpha_channel = min.(alpha.(img), masking_alpha_channel)
        return alphacolor.(img, combined_masking_alpha_channel)
    end
    return _apply
end

function masker(mask::AbstractArray, img::AbstractArray{<:Colorant})
    return masker(mask)(img)
end

function apply_mask(mask::AbstractArray, img::AbstractArray{<:Colorant})
    masking_alpha_channel = _mask_to_alpha(mask)
    if typeof(eltype(img)) <: TransparentColor
        # combine the mask with any existing mask on the image
        masking_alpha_channel = min.(alpha.(img), masking_alpha_channel)
    end
    return alphacolor.(img, masking_alpha_channel)
end

"""
    _mask_to_alpha(mask::AbstractArray)

Convert `mask` into an alpha channel to be applied to an image. 
"""
function _mask_to_alpha(mask::AbstractArray{<:TransparentColor})
    return ones(size(mask)) - alpha.(mask)
end

function _mask_to_alpha(mask::AbstractArray{<:AbstractRGB})
    return _mask_to_alpha(Gray.(mask))
end

function _mask_to_alpha(mask::AbstractArray{<:AbstractGray})
    return ones(size(mask)) - (gray.(mask) .> 0)
end

function _mask_to_alpha(mask::BitMatrix)
    return ones(size(mask)) - mask
end

function _mask_to_alpha(mask::AbstractArray{<:Real})
    return ones(size(mask)) - clamp.(mask, 0, 1)
end
