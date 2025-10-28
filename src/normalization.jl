
# TODO: Remove function
function imbinarize(img)
    f = AdaptiveThreshold(img) # infer the best `window_size` using `img`
    return binarize(img, f)
end
