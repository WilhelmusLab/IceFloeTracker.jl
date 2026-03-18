# Segmentation
The IFT segmentation functions include functions for semantic segmentation (pixel-by-pixel assignment into predefined categories) and object-based segmentation (groupings of pixels into distinct objects). Most of the segmentation algorithms currently focus on producing a binary mask, dividing objects of interest from the background. The semantic segmentation steps use $k$-means to group pixels into water and ice regions. A combination of watershed functions, morphological operations, and further applications of $k$-means are used for producing an object-based segmentation. 

## ``k``-means Binarization
``k``-means clustering is a common technique where points are grouped into a predetermined number of clusters. The IFT package contains a wrapper for the `Clustering.jl` ``k``-means algorithm enabling its use on grayscale images and arrays. Because sea ice floes tend to contain the brightest pixels in a scene, ``k``-means clustering is often effective at grouping ice floes together into a single cluster. If clouds and land have been masked, then the remaining pixels tend to come from only a few surface types -- water, new ice, and mixed sub-pixel-resolution ice floes. The ``k``-means binarization technique operates as follows:
1. Use ``k``-means clustering to identify coherent color regions in an image
2. Identify pixels with bright ice 
3. Select the cluster with the largest fraction of bright ice pixels.
The process is illustrated in the image below, as well as in an example notebook in the Github repository.

```@raw html
<img src="../assets/kmeans_example_case_006.png" width="600" alt="Example of k-means workflow. Shows a truecolor image of sea ice, a 4-color k-means segmentation, and two binarized images"/>
```

In the top left, we see the truecolor image for Case 006, from the Aqua satellite. Prior to the $k$-means clustering, we cast the image to grayscale, equalized it, and sharpened it. The top right shows the $k$-means clusters with $k=4$. In the bottom left is the result from the IceDetectionAlgorithm (here, `IceDetectionBrightnessPeaksMODIS721`). Finally, the $k$-means binarization result is in the bottom right.
