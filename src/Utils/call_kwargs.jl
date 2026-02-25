
"""
    call_kwargs(; name_to_function_map...)

Factory function for a callback which passes keyword arguments to functions specified in `name_to_function_map`. 
The callback can be used to pass intermediate results from a segmentation algorithm to user-specified functions, such as saving to disk or logging.
"""
function call_kwargs(; name_to_function_map...)
    return function callback(; kwargs...)
        for (name, func) in name_to_function_map
            func(kwargs[name])
        end
    end
    return callback
end
