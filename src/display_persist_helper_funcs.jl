# Helper functions
    """
        `make_filename()`

    Makes default filename with timestamp.
    
    """
    function make_filename()
        "persisted_mask-" * Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") * ".png"
    end


    function check_call(func_call::Expr)
        # Check expression is a function call or macrocall
        @assert func_call.head in [:call, :macrocall]  "Expected function call or macrocall; got a $(func_call.head)"
    end

    """
        `check_matrix(img)`

    Checks `img` is of type `Matrix`

    # Arguments
    - `img`: Image object (Matrix)
    """
    function check_matrix(img::Matrix)
        @assert img isa Matrix "Expected a Matrix type but got a $(typeof(img))."
    end

    """
        `check_fname(fname)`

    Checks `fname` does not exist in current directory; throws an assertion if this condition is false.

    # Arguments
    - `fname`: String object or Symbol to a reference to a String representing a path.
    """
    function check_fname(fname::Union{String, Symbol, Nothing})

        if fname isa String
            # println(:(fname)," is a type of ", typeof(fname), " and refers to the object ", fname)
            check_name = fname 
        elseif fname isa Symbol
            # println(:fname, " is a type of ", typeof(fname), " and refers to the object ", fname)
            check_name = eval(fname)
        elseif isnothing(fname)
            check_name = make_filename()
        end

        # println(fname, " is a ", typeof(fname))
        # println(:fname, " will turn into ", check_name)
        #     name_to_check = :($fname)
        # elseif typeof(fname) == Symbol
        #     println(fname, " is a ", typeof(fname))
        #     name_to_check = fname
        # end
        # return name_to_check
        @assert !isfile(check_name) "$check_name already exists in $(pwd())"
        return check_name
    end



    # # local tests
    # using Dates
    # bad_filename = "cameraman.png" # 
    # good_filename = "persisted_img-"*Dates.format(Dates.now(),"yyyy-mm-dd-HH:MM:SS:ss")*".png"
    # this_is_a_symbol = "pathto.file"
    # sym = :this_is_a_symbol
    # str = :"this is a string"
    # check_fname(sym)
    # check_fname(str)
    # check_fname(bad_filename)
    # check_fname(good_filename)
    # check_fname(nothing)