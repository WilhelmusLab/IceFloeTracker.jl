module skimage

using PyCall

export sk_morphology

const sk_morphology = PyNULL()

function __init__()
    skimage = "scikit-image=0.25.2"
    copy!(sk_morphology, pyimport_conda("skimage.morphology", skimage))
    return nothing
end

end
