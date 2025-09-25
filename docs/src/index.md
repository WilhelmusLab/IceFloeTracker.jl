# IceFloeTracker.jl

## Overview
IceFloeTracker.jl is a collection of routines and tools for processing remote sensing imagery, identifying sea ice floes, and tracking the displacement and rotation of ice floes across multiple images. It can be used either standalone to create custom processing pathways or with the [Ice Floe Tracker Pipeline](https://github.com/WilhelmusLab/ice-floe-tracker-pipeline).

```@contents
```

## Algorithm components
The Ice Floe Tracker (IFT) package includes the core functions for the three main steps of the algorithm. These functions can be used independently and can be customized for specific use cases. 

### Preprocessing
IFT operates on optical satellite imagery. The main functions are designed with "true color" and "false color" imagery in mind, and have thus far primarily been tested on imagery from the Moderate Resolution Imaging Spectroradiometer (MODIS) from the NASA _Aqua_ and _Terra_ satellites. The preprocessing routines mask land and cloud features, and aim to adjust and sharpen the remainder of the images to amplify the contrast along the edges of sea ice floes. (TBD: Link to main preprocessing page)

### Segmentation
The IFT segmentation functions include functions for semantic segmentation (pixel-by-pixel assignment into predefined categories) and object-based segmentation (groupings of pixels into distinct objects). The semantic segmentation steps use $k$-means to group pixels into water and ice regions. A combination of watershed functions, morphological operations, and further applications of $k$-means are used to identify candidate ice floes. (TBD: Link to main segmentation page)

### Tracking
Ice floe tracking is carried out by comparing the shapes produced in the segmentation step. Shapes with similar area are rotated until the difference in surface area is minimized, and then the edge shapes are compared using a Ѱ-s curve. If thresholds for correlation and area differences are met, then the floe with the best correlation and smallest area differences are considered matches and the objects are assigned the same label. In the end, trajectories for individual floes are recorded in a dataframe.



## Developers
IceFloeTracker.jl is a product of the [Wilhelmus Lab](https://www.wilhelmuslab.me) at Brown University, led by Monica M. Wilhelmus. The original algorithm was developed by Rosalinda Lopez-Acosta during her PhD work at University of California Riverside, advised by Dr. Wilhelmus. The translation of the original Matlab code into the current modular, open source Julia package has been carried out in conjunction with the Center for Computing and Visualization at Brown University. Contributors include Daniel Watkins, Maria Isabel Restrepo, Carlos Paniagua, Tim Divoll, John Holland, and Bradford Roarr.

## Citing

If you use IceFloeTracker.jl in research, teaching, or elsewhere, please mention the IceFloeTracker package and cite our journal article outlining the algorithm:

Lopez-Acosta et al., (2019). Ice Floe Tracker: An algorithm to automatically retrieve Lagrangian trajectories via feature matching from moderate-resolution visual imagery. _Remote Sensing of Environment_, **234(111406)**, doi:[10.1016/j.rse.2019.111406](https://doi.org/10.1016/j.rse.2019.111406).

## Papers using Ice Floe Tracker
1. Manucharyan, Lopez-Acosta, and Wilhelmus (2022)\*. Spinning ice floes reveal intensification of mesoscale eddies in the western Arctic Ocean. _Scientific Reports_, **12(7070)**, doi:[10.1038/s41598-022-10712-z](https://doi.org/10.1038/s41598-022-10712-z)
2. Watkins, Bliss, Hutchings, and Wilhelmus (2023)\*. Evidence of Abrupt Transitions Between Sea Ice Dynamical Regimes in the East Greenland Marginal Ice Zone. _Geophysical Research Letters_, **50(e2023GL103558)**, pp. 1-10, doi:[10.1029/2023GL103558](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023GL103558)

\*Papers using data from the Matlab implementation of Ice Floe Tracker.

## Functions

### Preprocessing

```@autodocs
Modules = [IceFloeTracker.Preprocessing]
Order   = [:function, :macro, :type]
Private = false
```

### Segmentation

```@autodocs
Modules = [IceFloeTracker.Segmentation]
Order   = [:function, :macro, :type]
Private = false
```

### Tracking

```@autodocs
Modules = [IceFloeTracker.Tracking]
Order   = [:function, :macro, :type]
Private = false
```

### Morphology

```@autodocs
Modules = [IceFloeTracker.Morphology]
Order   = [:function, :macro, :type]
Private = false
```

### Filtering

```@autodocs
Modules = [IceFloeTracker.Filtering]
Order   = [:function, :macro, :type]
Private = false
```

### Utils

```@autodocs
Modules = [IceFloeTracker.Utils]
Order   = [:function, :macro, :type]
Private = false
```

### Unsorted Functions

!!! todo "Functions here still need to be sorted"
    The functions which are shown in this section 
    will be reorganized into submodules.

```@autodocs
Modules = [IceFloeTracker]
Order   = [:module, :function, :macro, :type]
Private = false
```

## Index
```@index
```
