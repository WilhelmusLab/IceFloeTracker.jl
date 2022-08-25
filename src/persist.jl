# Persist macros
"""
    @persist img fname or
    @persist(img,fname) or
    @persist img or
    @persist(img)

Given a reference to an image object `img`, the macro persists (saves to a file) `img` to the current working directory using `fname` as filename. Returns `img`.

# Arguments
- `img`: Symbol expression representing an image object loaded in memory.
- `fname`: Optional filename for the persisted image.
"""
macro persist(_img, fname::Union{String,Symbol,Expr,Nothing}=nothing)
    return quote
        img = $(esc(_img))
        fname = check_fname($(esc(fname)))
        @info "Persisting image to file $(fname) in directory $(pwd()).\nTo load the persisted object use `load(img_path)`"
        save(fname, img)
        img
    end
end
