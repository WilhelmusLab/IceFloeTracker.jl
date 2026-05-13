"""
    call_kwargs(; name_to_function_map...)

Factory function for a callback which passes keyword arguments to functions specified in `name_to_function_map`. 
The callback can be used to pass intermediate results from a segmentation algorithm to user-specified functions, such as saving to disk or logging.

In this simple example, push_to_mutable_array is a simple call-back that pushes a value to a mutable array:

```julia-repl
julia> using IceFloeTracker.Utils: call_kwargs
julia> mutable = [] 
julia> function push_to_mutable_array(x)
           push!(mutable, x)
           return nothing
       end
julia> callback = call_kwargs(; value=push_to_mutable_array)
julia> callback(; value=5)
julia> mutable
1-element Vector{Any}:
 5
```

A more realistic example might have side-effects, 
such as printing a line, 
or saving intermediate segmentation results to disk:
```julia-repl
julia> using IceFloeTracker.Utils: call_kwargs
julia> function echo(x)
           println("Echoing: ", x)
           return nothing
       end
julia> callback = call_kwargs(; message=echo)
julia> callback(; message="[Segmented image data]")
Echoing: [Segmented image data]
```

"""
function call_kwargs(; name_to_function_map...)
    return function callback(; kwargs...)
        for (name, func) in name_to_function_map
            func(kwargs[name])
        end
    end
end
