include("display_persist_helper_funcs.jl")

# Main functions
"""
    `@display img`

Given an image object `img`, display `img` for visual inspection using `ImageView.imshow`

# Arguments
- `img`: Image object of type Matrix.
"""
macro display(img::Symbol)
    return quote

    # Display img
    output = $(esc(img))
    
    # Check output is a matrix
    check_matrix(output)

    println("Displaying image...")
    ImageView.imshow(output)
    output
    end
end

"""
    `@display create_mask_func(img)`

Given a function call `create_mask_func(img)` to build a mask (for functions such as `create_cloudmask` 
or `create_landmask`), the macro adds the following side effect to the function call:
- Display the mask for visual inspection using `ImageView.imshow`

# Arguments
- `create_mask_func(img)`: unevaluated function call expression with function `create_mask_func` and argument `img`.
"""
macro display(func_call::Expr)
    # Check expression is a function call or macrocall
    check_call(:($func_call))
    
    return quote 
        
        # Get function call output to return later 
        output = $(esc(func_call))
        
        # Display output
        println("Displaying image...")
        ImageView.imshow(output)

        output
    end
end

#  # Local test during development -- All good!
# using Images
# using Dates
# # using Pkg; Pkg.add("ImageView")
# using ImageView
# test_data_dir = "./test/data"
# img_path = "/landmask.tiff"
# img = test_data_dir * img_path |> Images.load
# @display img
# @persist img
# @persist @display img