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
