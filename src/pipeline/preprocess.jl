function check_landmask_path(lmpath::String)::Nothing
    name = basename(lmpath)
    input = dirname(lmpath)
    !isfile(lmpath) && error(
        "`$(name)` not found in $input. Please ensure a coastline image file named `$name` exists in $input.",
    )
    return nothing
end

"""
    cache_vector(type::Type, numel::Int64, size::Tuple{Int64, Int64})::Vector{type}

Build a vector of types `type` with `numel` elements of size `size`.

Example

```jldoctest
julia> cache_vector(Matrix{Float64}, 3, (2, 2))
3-element Vector{Matrix{Float64}}:
 [0.0 6.9525705991269e-310; 6.9525705991269e-310 0.0]
 [0.0 6.9525705991269e-310; 6.9525705991269e-310 0.0]
 [0.0 6.95257028858726e-310; 6.95257029000147e-310 0.0]
```
"""
function cache_vector(type::Type, numel::Int64, size::Tuple{Int64,Int64})::Vector{type}
    return [type(undef, size) for _ in 1:numel]
end

"""
    load_imgs(; input::String, image_type::String)

Load all images of type `image_type` (either `"truecolor"` or `"reflectance"`) in `input` into a vector.
"""
function load_imgs(; input::String, image_type::Union{Symbol,String})
    return [
        float64.(load(joinpath(input, f))) for
        f in readdir(input) if contains(f, string(image_type))
    ]
end

function load_truecolor_imgs(; input::String)
    return load_imgs(; input=input, image_type="truecolor")
end

function load_reflectance_imgs(; input::String)
    return load_imgs(; input=input, image_type="reflectance")
end

"""
    sharpen(truecolor_imgs::Vector{Matrix{Float64}}, landmask_no_dilate::Matrix{Bool})

Sharpen truecolor images with the non-dilated landmask applied. Returns a vector of sharpened images.

"""
function sharpen(
    truecolor_imgs::Vector{Matrix{RGB{Float64}}}, landmask_no_dilate::BitMatrix
)
    @info "Sharpening truecolor images..."
    return [imsharpen(img, landmask_no_dilate) for img in truecolor_imgs]
end

function cloudmask(; input::String, output::String)::Vector{BitMatrix}
    # find reflectance imgs in input dir
    ref = [img for img in readdir(input) if contains(img, "reflectance")] # ref is sorted
    total_ref = length(ref)
    @info "Found $(total_ref) reflectance images in $input. 
    Cloudmasking false color images..."

    # Preallocate container for the cloudmasks
    ref_img = IceFloeTracker.float64.(IceFloeTracker.load(joinpath(input, ref[1]))) # read in the first one to retrieve size
    sz = size(ref_img)
    cloudmasks = [BitMatrix(undef, sz) for _ in 1:total_ref]

    # Do the first one because it's loaded already
    cloudmasks[1] = IceFloeTracker.create_cloudmask(ref_img)
    # and now the rest
    for i in 2:total_ref
        img = IceFloeTracker.float64.(IceFloeTracker.load(joinpath(input, ref[i])))
        cloudmasks[i] = IceFloeTracker.create_cloudmask(img)
    end
    return cloudmasks
end

"""
    disc_ice_water(
    reflectance_imgs::Vector{Matrix{RGB{Float64}}},
    normalized_imgs::Vector{Matrix{Gray{Float64}}},
    cloudmasks::Vector{BitMatrix},
    landmask::BitMatrix,
)

Generate vector of ice/water discriminated images from the collection of reflectance, sharpened truecolor, and cloudmask images using the study area landmask. Returns a vector of ice/water masks.

"""
function disc_ice_water(
    reflectance_imgs::Vector{Matrix{RGB{Float64}}},
    normalized_imgs::Vector{Matrix{Gray{Float64}}},
    cloudmasks::Vector{BitMatrix},
    landmask::BitMatrix,
)
    return [
        IceFloeTracker.discriminate_ice_water(ref_img, norm_img, landmask, cldmsk) for
        (ref_img, norm_img, cldmsk) in zip(reflectance_imgs, normalized_imgs, cloudmasks)
    ]
end

"""
    sharpen_gray(
    sharpened_imgs::Vector{Matrix{Float64}},
    landmask::AbstractArray{Bool},
)

Apply the landmask to the collection of sharpened truecolor images and return a gray colorview of the collection.
"""
function sharpen_gray(
    sharpened_imgs::Vector{Matrix{Float64}},
    landmask::AbstractArray{Bool},
)
    return [IceFloeTracker.imsharpen_gray(img, landmask) for img in sharpened_imgs]
end

function get_ice_labels(
    reflectance_imgs::Vector{Matrix{RGB{Float64}}},
    landmask::AbstractArray{Bool}
)
    return [IceFloeTracker.find_ice_labels(ref_img, landmask) for ref_img in reflectance_imgs]
end

"""
    preprocess(; truecolor_image, reflectance_image, landmask_imgs)

Preprocess and segment floes in `truecolor_image` and `reflectance_image` images using the landmasks  `landmask_imgs`. Returns a boolean matrix with segmented floes for feature extraction.

# Arguments
- `truecolor_image::T`: truecolor image to be processed
- `reflectance_image::T`: reflectance image to be processed
- `landmask_imgs`: named tuple with dilated and non-dilated landmask images
"""
function preprocess(truecolor_image::T, reflectance_image::T, landmask_imgs::NamedTuple{(:dilated, :non_dilated),Tuple{BitMatrix,BitMatrix}}) where {T<:Matrix{RGB{Float64}}}

    @info "Building cloudmask"
    cloudmask = create_cloudmask(reflectance_image)

    # 2. Intermediate images
    @info "Finding ice labels"
    ice_labels = IceFloeTracker.find_ice_labels(reflectance_image, landmask_imgs.non_dilated)

    @info "Sharpening truecolor image"
    # a. apply imsharpen to truecolor image using non-dilated landmask
    sharpened_truecolor_image = IceFloeTracker.imsharpen(truecolor_image, landmask_imgs.non_dilated)
    # b. apply imsharpen to sharpened truecolor img using dilated landmask
    sharpened_gray_truecolor_image = IceFloeTracker.imsharpen_gray(sharpened_truecolor_image, landmask_imgs.dilated)

    @info "Normalizing truecolor image"
    normalized_image = IceFloeTracker.normalize_image(
        sharpened_truecolor_image, sharpened_gray_truecolor_image, landmask_imgs.dilated)

    # Discriminate ice/water
    @info "Discriminating ice/water"
    ice_water_discrim = IceFloeTracker.discriminate_ice_water(
        reflectance_image, normalized_image, copy(landmask_imgs.dilated), cloudmask)

    # 3. Segmentation
    @info "Segmenting floes part 1/3"
    segA = IceFloeTracker.segmentation_A(IceFloeTracker.segmented_ice_cloudmasking(
        ice_water_discrim, cloudmask, ice_labels
    ))

    # segmentation_B
    @info "Segmenting floes part 2/3"
    segB = IceFloeTracker.segmentation_B(sharpened_gray_truecolor_image, cloudmask, segA)

    # Process watershed in parallel using Folds
    @info "Building watersheds"

    watersheds_segB = pmap(IceFloeTracker.watershed_ice_floes, [segB.not_ice_bit, segB.ice_intersect])
    # reuse the memory allocated for the first watershed
    watersheds_segB[1] .= IceFloeTracker.watershed_product(watersheds_segB...)

    # segmentation_F
    @info "Segmenting floes 3/3"
    return IceFloeTracker.segmentation_F(
        segB.not_ice,
        segB.ice_intersect,
        watersheds_segB[1],
        ice_labels,
        cloudmask,
        landmask_imgs.dilated,
    )
end

"""
    preprocess(truedir::T, refdir::T, lmdir::Tuple{T, T}, output::T)

Preprocess and segment floes in all images in `truedir` and `refdir` using the landmasks in `lmdir`. Saves the segmented floes in `output`.

# Arguments
- `truedir`: directory with truecolor images to be processed
- `refdir`: directory with reflectance images to be processed
- `lmdir`: directory with dilated and non-dilated landmask images
"""
function preprocess(; truedir::T, refdir::T, lmdir::T, output::T) where {T<:AbstractString}
    # 1. Get references to images
    truecolor_refs = [ref for ref in readdir(truedir) if occursin("truecolor", ref)]
    reflectance_refs = [ref for ref in readdir(refdir) if occursin("reflectance", ref)]
    landmask_imgs = deserialize(joinpath(lmdir, "generated_landmask.jls"))

    # 2. Preprocess
    @info "Preprocessing"
    truecolor_container = loadimg(dir=truedir, fname=truecolor_refs[1]) |> similar
    reflectance_container = similar(truecolor_container)
    segmented_floes = cache_vector(BitMatrix, length(truecolor_refs), size(truecolor_container))
    numimgs = length(truecolor_refs)
    Threads.@threads for i in eachindex(truecolor_refs)
        @info "Processing image $i of $numimgs"
        truecolor_container .= loadimg(dir=truedir, fname=truecolor_refs[i])
        reflectance_container .= loadimg(dir=refdir, fname=reflectance_refs[i])
        segmented_floes[i] .= preprocess(truecolor_container, reflectance_container, landmask_imgs)
    end

    # 3. Save
    @info "Serializing segmented floes"
    serialize(joinpath(output, "segmented_floes.jld"), segmented_floes)
    nothing
end
