function check_landmask_path(lmpath::String)::Nothing
    name = basename(lmpath)
    input = dirname(lmpath)
    !isfile(lmpath) && error(
        "`$(name)` not found in $input. Please ensure a coastline image file named `$name` exists in $input.",
    )
    return nothing
end

"""
    landmask(; input, output)

Given an input directory with a landmask file and possibly truecolor images, create a land/soft ice mask. The resulting images are saved to the snakemake output directory. Returns the landmkask object. 

# Arguments
- `input`: path to image dir containing truecolor and landmask source images
- `output`: path to output dir where land-masked truecolor images and the generated binary land mask are saved

"""
function landmask(; input::String, output::String)
    landmask_fname = "landmask.tiff"
    @info "Looking for $landmask_fname in $input"

    lmpath = joinpath(input, landmask_fname)
    check_landmask_path(lmpath)
    @info "$landmask_fname found in $input. Creating landmask..."

    img = load(lmpath)
    mkpath(output)
    out = @persist create_landmask(img) joinpath(output, "generated_landmask.png")
    @info "Landmask created succesfully."
    return out
end

function cloudmask_reflectance(input::String, output::String)::Vector{BitMatrix}
    # find reflectance imgs in input dir
    ref = sort([img for img in readdir(input) if contains(img, "reflectance")])
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
        cloudmasks[i] = IceFloeTracker.create_cloudmask(ref_img)
    end
    return cloudmasks
end
