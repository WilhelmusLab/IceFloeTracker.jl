"""
    apply_to_channels(img, f)

Broadcasts a function f to each channel of input image img, then recombines to return.

"""
function apply_to_channels(img, f)
    _view = channelview(img)
    _result = [f(@view(_view[i, :, :])) for i in 1:3]
    return colorview(eltype(img), _result...) 
end
