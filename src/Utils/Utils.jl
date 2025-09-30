module Utils

export @persist

using Images: save, Gray, FixedPoint
using Dates: format, now

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
        save(fname, img)
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

"""
    check_fname(fname)

Checks `fname` does not exist in current directory; throws an assertion if this condition is false.

# Arguments
- `fname`: String object or Symbol to a reference to a String representing a path.
"""
function check_fname(fname::Union{String,Symbol,Nothing}=nothing)::String
    if fname isa String # then use as filename
        check_name = fname
    elseif fname isa Symbol
        check_name = eval(fname) # get the object represented by the symbol
    elseif isnothing(fname) # nothing provided so make a filename
        check_name = make_filename()
    end

    # check name does not exist in wd
    isfile(check_name) && error("$check_name already exists in $(pwd())")
    return check_name
end

"""
    make_filename()

Makes default filename with timestamp.

"""
function make_filename()::String
    return "persisted_img-" * format(now(), "yyyy-mm-dd-HHMMSS") * ".png"
end

"""
    timestamp(fname)

Attach timestamp to `fname`.
"""
function timestamp(fname::String)
    ts = format(now(), "yyyy-mm-dd-HHMMSS")
    return fname * "-" * ts
end

"""
    fname_ext_split(fname)

Split `"fname.ext"` into `"fname"` and `"ext"`.
"""
function fname_ext_split(fname::String)
    return (name=fname[1:(end - 4)], ext=fname[(end - 2):end])
end

"""
    fname_ext_splice(fname, ext)

Join `"fname"` and `"ext"` with `'.'`.
"""
function fname_ext_splice(fname::String, ext::String)
    return fname * '.' * ext
end

end
