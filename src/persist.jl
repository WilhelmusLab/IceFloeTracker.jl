# Persist macros
"""
    @persist img fname
    @persist(img,fname)
    @persist img
    @persist(img)
    @persist img fname ts
    @persist(img, fname, ts)

Given a reference to an image object `img`, the macro persists (saves to a file) `img` to the current working directory using `fname` as filename. Returns `img`.

# Arguments
- `img`: Symbol expression representing an image object loaded in memory.
- `fname`: Optional filename for the persisted image.
- `ts`: Optional boolean to attach timestamp to `fname`.
"""
macro persist(_img, fname::Union{String,Symbol,Expr,Nothing}=nothing)
    return quote
        img = $(esc(_img))
        fname = check_fname($(esc(fname)))
        @info "Persisting image to $(fname).\nTo load the persisted object use `Images.load(img_path)`"
        Images.save(fname, img)
        img
    end
end

macro persist(_img, _fname::Union{String,Symbol,Expr,Nothing}, ts::Bool)
    # fname provided and ts desired
    if !isnothing(_fname) && ts
        _fname = string(_fname)
        if _fname[end - 3] == '.' # with an ext?
            name, ext = fname_ext_split(_fname)
            # tack on the timestamp and splice the ext back
            _fname = fname_ext_splice(timestamp(name), ext)
        end
    end
    return quote
        @persist($(esc(_img)), $(esc(_fname)))
    end
end
