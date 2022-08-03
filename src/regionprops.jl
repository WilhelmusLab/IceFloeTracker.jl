# remove these later
using PyCall
using DataFrames

# PyCall imports
skimage = pyimport("skimage")
sk_regionprops = skimage.measure.regionprops
sk_label = skimage.measure.label
sk_regionprops_table = skimage.measure.regionprops_table

"""
    regionprops(bw_img)

A wrapper of the `regionprops_table` function from the skimage python library. See its documentation at https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.regionprops.
    

"""
function regionprops(bw_img::Any;
                     properties::Vector{String},
                     connectivity::Int=2)

    # bw_img read in with Images.load
    if eltype(bw_img) != Float64
        bw_img = Float64.(bw_img)
    end

    labels = sk_label(bw_img, connectivity)

    
    return DataFrame(sk_regionprops_table(labels, bw_img, properties))


end
    

"""
    to_dataframe(dict)

Convert a Dict type object (such as the output of `regionprops_table`) to a DataFrame (table-like) type object.
"""
function to_dataframe(dict::Dict)
    return DataFrame(dict)
    
end
"""
    sk_label(img)
"""
function sk_label(img, background=None, return_num=False, connectivity=None)
    
end
