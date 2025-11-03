module skimage

using PythonCall

export sk_measure, sk_morphology, sk_exposure

const sk_measure = PythonCall.pynew()
const sk_morphology = PythonCall.pynew()
const sk_exposure = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(sk_measure, pyimport("skimage.measure"))
    PythonCall.pycopy!(sk_morphology, pyimport("skimage.morphology"))
    PythonCall.pycopy!(sk_exposure, pyimport("skimage.exposure"))
end

end
