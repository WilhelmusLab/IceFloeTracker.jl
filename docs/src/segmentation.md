# Segmentation
The IFT segmentation functions include functions for semantic segmentation (pixel-by-pixel assignment into predefined categories) and object-based segmentation (groupings of pixels into distinct objects). Most of the segmentation algorithms currently focus on producing a binary mask, dividing objects of interest from the background. The semantic segmentation steps use $k$-means to group pixels into water and ice regions. A combination of watershed functions, morphological operations, and further applications of $k$-means are used for producing an object-based segmentation. 

## ``k``-means Binarization
``k``-means clustering is a common technique where points are grouped into a predetermined number of clusters. The IFT package contains a wrapper for the `Clustering.jl` ``k``-means algorithm enabling its use on grayscale images and arrays. Because sea ice floes tend to contain the brightest pixels in a scene, ``k``-means clustering is often effective at grouping ice floes together into a single cluster. If clouds and land have been masked, then the remaining pixels tend to come from only a few surface types -- water, new ice, and mixed sub-pixel-resolution ice floes. The ``k``-means binarization technique operates as follows:
1. Use ``k``-means clustering to identify coherent color regions in an image
2. Identify pixels with bright ice 
3. Select the cluster with the largest fraction of bright ice pixels.
The process is illustrated in the image below, as well as in an example notebook in the Github repository.

