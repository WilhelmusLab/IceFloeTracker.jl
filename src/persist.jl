include("display_persist_helper_funcs.jl")

# Persist macros
"""
    `@persist img fname` or
    `@persist(img,fname)` or
    `@persist img` or
    `@persist(img)`

Given <`img::Symbol`> refering to an image object `img`, the macro persists (saves to a file) `img` to the current working directory using <`fname`> as filename.

# Arguments
- `img`: Symbol expression representing an image object loaded in memory.
- `fname`: Optional filename for the persisted image.
"""
macro persist(img::Symbol,
              fname::Union{String,Symbol,Nothing}=nothing)
    return quote
        # check img is an image object (Matrix)
        check_matrix($(esc(img)))
        # local fname = check_fname($(esc(fname)))
        fname = check_fname($(esc(fname)))
        local msg = "Persisting image to file $(fname) in directory $(pwd())"
        # msg = "Persisting image to file $($fname) in directory $(pwd())"
        println(msg)
        # @info msg
        # @info $msg
    #     # println("To load the persisted object use `JLD2.load_object(object)`")
        println("To load the persisted object use `load(img_path)`")
    #     # JLD2.save_object(filename, output)
        # Images.save($fname, $(esc(img)))
        Images.save(fname, $(esc(img)))
    #     println("Object persisted successfully to\n",filename)
    
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