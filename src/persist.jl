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
        @info "Persisting image to file $(fname) in directory $(pwd()).\nTo load the persisted object use `Images.load(img_path)`"
        Images.save(fname, img)
        img
    end
end

macro persist(_img, fname::Union{String,Symbol,Expr,Nothing}, ts::Bool=false)
    return quote
        img = $(esc(_img))
        fname = check_fname($(esc(fname)))
        ts && (fname = timestamp(fname))
        @info "Persisting image to file $(fname) in directory $(pwd()).\nTo load the persisted object use `Images.load(img_path)`"
        Images.save(fname, img)
        img
    end
end

function timestamp(fname::String)::String
    return fname * "_" * Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS")
end