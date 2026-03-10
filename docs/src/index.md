# IceFloeTracker.jl

## Overview
IceFloeTracker.jl is a collection of routines and tools for processing remote sensing imagery, identifying sea ice floes, and tracking the displacement and rotation of ice floes across multiple images. It can be used either standalone to create custom processing pathways or with the [Ice Floe Tracker Pipeline](https://github.com/WilhelmusLab/ice-floe-tracker-pipeline).

## Developers
IceFloeTracker.jl is a product of the [Wilhelmus Lab](https://www.wilhelmuslab.me) at Brown University, led by Monica M. Wilhelmus. The original algorithm was developed by Rosalinda Lopez-Acosta during her PhD work at University of California Riverside, advised by Dr. Wilhelmus. The translation of the original Matlab code into the current modular, open source Julia package has been carried out in conjunction with the Center for Computing and Visualization at Brown University. Contributors include Daniel Watkins, Maria Isabel Restrepo, Carlos Paniagua, Tim Divoll, John Holland, and Bradford Roarr.

## Citing
If you use IceFloeTracker.jl in research, teaching, or elsewhere, please mention the IceFloeTracker package and cite our journal article outlining the algorithm:

Lopez-Acosta et al., (2019). Ice Floe Tracker: An algorithm to automatically retrieve Lagrangian trajectories via feature matching from moderate-resolution visual imagery. _Remote Sensing of Environment_, **234(111406)**, doi:[10.1016/j.rse.2019.111406](https://doi.org/10.1016/j.rse.2019.111406).

## Papers using Ice Floe Tracker
1. Manucharyan, Lopez-Acosta, and Wilhelmus (2022)\*. Spinning ice floes reveal intensification of mesoscale eddies in the western Arctic Ocean. _Scientific Reports_, **12(7070)**, doi:[10.1038/s41598-022-10712-z](https://doi.org/10.1038/s41598-022-10712-z)
2. Covington, Chen, and Wilhelmus (2022)\*. Bridging Gaps in the Climate Observation Network: A Physics‐based Nonlinear Dynamical Interpolation of Lagrangian Ice Floe Measurements via Data‐Driven Stochastic Models. _Journal of Advances in Modeling Earth Systems_, **14 (e2022MS003218)**, pp. 1-28, doi:[10.1029/2022MS003218](https://doi.org/10.1029/2022MS003218)
3. Watkins, Bliss, Hutchings, and Wilhelmus (2023)\*. Evidence of Abrupt Transitions Between Sea Ice Dynamical Regimes in the East Greenland Marginal Ice Zone. _Geophysical Research Letters_, **50(e2023GL103558)**, pp. 1-10, doi:[10.1029/2023GL103558](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023GL103558)
4. Watkins, Buckley, Kim, Hutchings, and Wilhelmus (2025). Characterizing the Marginal Ice Zone in the Greenland Sea Through Seasonal Floe-Scale Sea Ice Observations. _ESS Open Archive (submitted to Journal of Glaciology)_. doi:[10.22541/essoar.175503380.04848000/v1](https://doi.org/10.22541/essoar.175503380.04848000/v1)

\*Papers using data from the Matlab implementation of Ice Floe Tracker.
