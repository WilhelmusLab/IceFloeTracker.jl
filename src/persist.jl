include("display_persist_helper_funcs.jl")

# Persist macros
"""
    `@persist img fname` or
    `@persist(img,fname)` or
    `@persist img` or
    `@persist(img)`

Given `img::Symbol` refering to an image object `img`, the macro persists (saves to a file) `img` to the current working directory using `fname` as filename. Returns the persisted image.

# Arguments
- `img`: Symbol expression representing an image object loaded in memory.
- `fname`: Optional filename for the persisted image.
"""
macro persist(_img::Symbol,
              fname::Union{String,Symbol,Nothing}=nothing)
    return quote
        img = $(esc(_img))
        # check img is an image object (Matrix)
        check_matrix(img)
        # local fname = check_fname($(esc(fname)))
        fname = check_fname($(esc(fname)))
        local msg = "Persisting image to file $(fname) in directory $(pwd())"
        # msg = "Persisting image to file $($fname) in directory $(pwd())"
        println(msg)
        
        println("To load the persisted object use `Images.load(img_path)`")
    
        Images.save(fname, img)
        println("Object persisted successfully to\n",fname)
        img
    end
end

"""
    `@persist create_mask_func(img) fname` or
    `@persist(create_mask_func(img), fname=nothing)` or
    `@persist create_mask_func(img)` or
    `@persist(create_mask_func(img))`


Given a function call `create_mask_func(img)` and an optional `filename` to build a mask (for functions such as `create_cloudmask` or `create_landmask`), the macro adds the following side effect to the function call:
- Persists the generated mask to an image file using `Images.save`.

# Arguments
- `create_mask_func(img)`: unevaluated function call expression with function `create_mask_func` and argument `img`.
- `fname`: Optional filename for the persisted image.
"""
macro persist(func_call::Expr,
              fname::Union{String,Symbol,Nothing}=nothing)
    
    # Check expression is a call
    check_call(:($func_call))

    return quote
    
    # # Get function call output to return later 
    img = $(esc(func_call))
    
    # # Persist output
    fname = check_fname($(esc(fname)))
    println("Persisting mask to file $fname in directory $(pwd())")
    println("To load the persisted object use `Images.load(object)`")
    Images.save(fname, img)
    println("Object persisted successfully to\n",fname)
    img
    end
end

# # Local test during development -- All good!
# using Dates
# using Images
# test_data_dir = "./test/data"
# img_path = "/landmask.tiff"
# outimage_path = "outimage1.tiff"
# img = test_data_dir * img_path |> Images.load
# # img = TestImages.testimage("camera");
# @persist img "image595i.png"
# @persist img outimage_path
# @persist img
# @assert isfile("image595i.png")
# @assert isfile(outimage_path)
# rm("image595i.png"); rm(outimage_path);
# for f in readdir()
#     if startswith(f,"persisted_mask")
#         rm(f)
#     end
# end
# # all good!

# # # # Local test during development 2
# using Dates
# using Images
# test_data_dir = "./test/data"
# img_path = "/landmask.tiff"
# outimage_path = "outimage1.tiff"
# img = test_data_dir * img_path |> Images.load
# # img = TestImages.testimage("camera");
# # @persist identity(img) "image595i.png"
# # @persist identity(img) outimage_path
# @persist identity(img)
# @assert isfile("image595i.png")
# @assert isfile(outimage_path)
# rm("image595i.png"); rm(outimage_path);
# for f in readdir()
#     if startswith(f,"persisted_mask")
#         rm(f)
#     end
# end
# all good!