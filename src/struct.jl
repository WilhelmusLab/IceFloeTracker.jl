"""
Turn a struct into a collection of ::Symbol => value pairs.
"""
_unpack_struct(x::Any) = [k => getproperty(x, k) for k in fieldnames(typeof(x))]
