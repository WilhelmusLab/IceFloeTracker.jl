# Helper functions
"""
    `make_filename()`

Makes default filename with timestamp.

"""
function make_filename()
    return "persisted_mask-" * Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") * ".png"
end

"""
    `check_fname(fname)`

Checks `fname` does not exist in current directory; throws an assertion if this condition is false.

# Arguments
- `fname`: String object or Symbol to a reference to a String representing a path.
"""
function check_fname(fname::Union{String,Symbol,Nothing}=nothing)
    if fname isa String # then use as filename
        check_name = fname
    elseif fname isa Symbol
        check_name = eval(fname) # get the object represented by the symbol
    elseif isnothing(fname) # nothing provided so make a filename
        check_name = make_filename()
    end

    # check name does not exist in wd
    @assert !isfile(check_name) "$check_name already exists in $(pwd())"
    return check_name
end

"""
    `add_padding(img, rad, type, val)`

Extrapolate the image `img` with `val` `rad` units beyond its boundary. Returns the extrapolated image.

# Arguments
- `img`: Image to be padded.
- `rad`: Uniform number of rows/columns to extrapolate beyond the image boundary.
- `type`: Symbol representing the type of extrapolation; defaults to `:replicate`. The other supported type is `type=:replicate`.
- `val`: Value to be used for the extrapolation (when `type=:fill`).
"""
function add_padding(img, rad::Int=0, type::Symbol=:replicate, val::Int=0)::Matrix
    if type == :replicate
        return collect(Images.padarray(img, Pad(type,rad,rad)))
    elseif type == :fill
        return collect(Images.padarray(img, Fill(val, (rad,rad),(rad,rad)))) 
    end
end

"""
    `remove_padding(img, rad)`

Removes `rad` units of padding uniformly along all sides of the image `img`. Returns the cropped image.

# Arguments
- `img`: Image to be padded.
- `rad`: Number of rows/columns to extrapolate beyond the image boundary. 
"""
function remove_padding(img, rad::Int)::Matrix
    return img[rad+1:end-rad,rad+1:end-rad]
end
