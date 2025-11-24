"""
Module for loading validated ice floe data.
"""
module Data

export AbstractLoader,
    AbstractLoader,
    GitHubLoader,
    Case,
    Dataset,
    loader,
    metadata,
    Watkins2026Dataset,
    metadata,
    name,
    modis_truecolor,
    modis_falsecolor,
    modis_landmask,
    modis_cloudfraction,
    masie_landmask,
    masie_seaice,
    validated_binary_floes,
    validated_labeled_floes,
    validated_floe_properties

include("./loader.jl")
include("./dataset.jl")
include("./watkins-2026.jl")

end # module Data
