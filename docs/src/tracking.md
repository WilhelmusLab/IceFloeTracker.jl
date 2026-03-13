# Tracking
Floe tracking in Ice Floe Tracker is based on object shapes and locations. Given a series of segmented images, the goal is to identify shapes that persist across images, subject to constraints on object similarity and maximum displacements. The algorithm has three main components: 
1. the main `FloeTracker` functor and internal `floe_tracker` function, which keeps track of the most recent floe in each trajectory (trajectory heads) and the set of candidate matches in the next image, 
2. a `FilterFunction` which identifies all the possible pairs between a trajectory head and floes in the candidate set,
3. a `MatchingFunction` which resolves conflicts to find pairs.

By designing and customizing the filter functions and matching functions, `FloeTracker` provides a flexible and powerful platform for developing floe tracking workflows.

## Preparing data for tracking
The current implementation of `FloeTracker` accepts two arguments: a vector of images and a vector of DateTimes. The vector of images may either be SegmentedImages or labeled integer arrays (background=0). The two vectors must be the same length, and each image must be the same size. For example, if `img_paths` contained an ordered list of saved, segmented binary images: 

```julia
using IceFloeTracker
using Images

binary_images = load.(img_paths)
labeled_images = label_components.(binary_images)

```

## Setting up the FloeTracker
The FloeTracker functor is initialized with a FilterFunction and a Matching function. Initializing with the default functions:
```julia
tracker = FloeTracker(FilterFunction(), MatchingFunction())
```

Optional keyword arguments specify the minimum and maximum floe sizes (in pixels) and the maximum time step in between floe pairs.
```julia
using Dates
minimum_area = 100
maximum_area = 90e3
maximum_time_step = Day(2)
tracker = FloeTracker(FilterFunction(), MatchingFunction(), minimum_area, maximum_area, maximum_time_step)
```

With image times `passtimes`, you can then track floes simply by running
```julia
tracked_floes = tracker(labeled_images, passtimes)
```

## Filter Functions
Floe tracking begins by feature extraction. `FloeTracker` uses the `regionprops_table` function to extract geometric information from the segmented
images. This function provides the following measures:
* label
* centroid (row, col)
* bounding box
* area
* convex area
* major and minor axis length
* perimeter
* orientation
In addition to the region props quantities, each label is given a random 12-character string UUID, a cropped binary mask with the segment shape, a computed ψ-s curve, and the image time. 

The `floe_tracker` takes three positional arguments: a list of DataFrames, a filter function, and a matching function, as  Assuming that a set of images `segmented_images` has already been produced, the `regionprops_table` function produces a DataFrame where each row corresponds to a floe, and each column is some measurement or attribute of the floe. DataFrames are highly flexible--entries in the dataframe are not limited to words and numbers, but can include vectors and matrices as well. At the very least, the property tables will need to include the floe ID, floe area and an associated observation time. Other columns and measures depend on what will be used in the filter function. The `regionprops_table` function defaults to calculating  area, convex area, centroid, perimeter, major and minor axis, and orientation. We can initialize the property tables using dot notation to broadcast to a list:

```julia
props = regionprops_table.(segmented_images)
```

We include helper functions to add unique IDs to each row and to add image observation times. Assuming `passtimes` is a list of DateTimes of the same length as `segmented_images`, we run

```julia
add_uuids!.(props)
add_passtimes!.(props, passtimes)
```

The default filter functions include calculations based on binary floe shapes (floe masks) and associated $\psi$-s curves. We include helper functions for these as well. For the floe masks, we also need a list with binary images associated with each segmented image.
```julia
add_floemasks!.(props, binary_images)
add_ψs!.(props)
```

## Floe Filter Functions
Floe filter functions take two argmuments: a DataFrameRow corresponding to the floe to be matched, and a DataFrame with candidate pairs from the current time step. These functions should operate in-place and result in DataFrame subsetting operation. IceFloeTracker.jl includes four main `AbstractFloeFilterFunctions`.  These functions all follow the same procedure:
1. Compute comparisons between the floe and each floe in the candidates DataFrame
2. Use an threshold test to evaluate the comparisons
3. Reduce the candidates dataframe to only those pairs that pass the threshold test
Each of these functions can be called in series, as they only depend on the columns in the input property tables. The`ChainedFilterFunction`takes a list of filter functions and wraps them into a single function call. Using the struct/functor approach, we can initialize each function with parameters, then pass the function and settings into the `floe_tracker` function. As an example, let's define a filter function which compares the relative error in area against a step function. We initialize the step function first:

```julia
sw = StepwiseLinearThresholdFunction(changepoint_area=500, low_value=0.25, high_value=0.125)
```
The function `sw` will return `true` if the relative error is less than 0.25 and the area is smaller than 500 pixels. If the area is larger, then a stricter threshold of 0.125 is applied. Next, we make a filter function
```julia
rel_err_ff = RelativeErrorThresholdFilter(variable=:area, threshold_function=sw)
```

This function will add a column `:relative_error_area` to the filter function, which can be used in the Matching Function.

## Matching Function


