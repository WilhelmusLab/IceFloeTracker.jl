module skimage

using PyCall

export sk_morphology, sk_exposure

const sk_morphology = PyNULL()
const sk_exposure = PyNULL()

function __init__()
    skimage = "scikit-image=0.25.2"
    copy!(sk_exposure, pyimport_conda("skimage.exposure", skimage))
    copy!(sk_morphology, pyimport_conda("skimage.morphology", skimage))
    return nothing
end

end
