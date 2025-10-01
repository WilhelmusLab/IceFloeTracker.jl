"""
Gradient Functions for Nonlinear Diffusion

This module defines the supported gradient functions for nonlinear diffusion algorithms,
particularly for the Perona-Malik diffusion implementation.
"""

# Define gradient functions
exponential(norm∇I, k) = exp(-(norm∇I / k)^2)
inverse_quadratic(norm∇I, k) = 1 / (1 + (norm∇I / k)^2)

struct SupportedFunctions
    functions::Dict{String,Function}

    function SupportedFunctions()
        funcs = Dict{String,Function}(
            "exponential" => exponential,
            "inverse_quadratic" => inverse_quadratic,
        )
        return new(funcs)
    end
end

# Allow it to behave like a collection for `in` operations
Base.in(x, sf::SupportedFunctions) = x in keys(sf.functions)

# Get the function implementation
Base.getindex(sf::SupportedFunctions, key::String) = sf.functions[key]

# Custom string representation for error messages
function Base.show(io::IO, sf::SupportedFunctions)
    return print(io, join(keys(sf.functions), ", "))
end

# Get all available function names
Base.keys(sf::SupportedFunctions) = keys(sf.functions)

# Check if a function is supported
function is_supported(sf::SupportedFunctions, name::String)
    return name in sf
end

const SUPPORTED_GRADIENT_FUNCTIONS = SupportedFunctions()
