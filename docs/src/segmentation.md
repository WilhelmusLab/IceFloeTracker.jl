# Segmentation
The IFT segmentation functions include functions for semantic segmentation (pixel-by-pixel assignment into predefined categories) and object-based segmentation (groupings of pixels into distinct objects). The semantic segmentation steps use $k$-means to group pixels into water and ice regions. A combination of watershed functions, morphological operations, and further applications of $k$-means are used to identify candidate ice floes. 

## Ice/Water Discrimination

## Feature Identification