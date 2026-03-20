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
In addition to the region props quantities, each label is given a random 12-character string UUID, a cropped binary mask with the segment shape, a computed ψ-s curve, and the image time. Each floe's information is stored in a row of a DataFrame.

## Filter Functions

Starting from the first floe ``f`` in the first image and corresponding time ``t``, `FloeTracker` selects all the DataFrame rows with time stamps up to `maximum_time_step` before and after ``t`` to form a set of candidate matches ``C``. The goal of the filter function is to reduce the size of ``C`` as far as possible without removing the true match (if it exists). 

Floe filter functions have type `AbstractFloeFilterFunction` and take two argmuments: a DataFrameRow corresponding to the floe to be matched, and a DataFrame with candidate pairs from the current time step. These functions operate in-place, modifying and subsetting DataFrame ``C``. IceFloeTracker.jl includes four main `AbstractFloeFilterFunctions`.  These functions all follow the same procedure:
1. Compute comparisons between the floe and each floe in the candidates DataFrame
2. Use an threshold test to evaluate the comparisons
3. Reduce the candidates dataframe to only those pairs that pass the threshold test
Each of these functions can be called in series, as they only depend on the columns in the input property tables. 

As an example, let's define a filter function which compares the relative error in area against a step function. We initialize the step function first:

```julia
sw = StepwiseLinearThresholdFunction(changepoint_area=500, low_value=0.25, high_value=0.125)
```
The function `sw` will return `true` if the relative error is less than 0.25 and the area is smaller than 500 pixels. If the area is larger, then a stricter threshold of 0.125 is applied. Next, we make a filter function
```julia
rel_err_ff = RelativeErrorThresholdFilter(variable=:area, threshold_function=sw)
```

This function will add a column `:relative_error_area` to the filter function, which can be used in the Matching Function.

Since the input and output of filter functions are dataframes, they can be linked together. We recommend considering the computational cost of each comparison metric, and ordering the tests from least expensive to most expensive. 

The `ChainedFilterFunction` functor takes a list of filter functions and wraps them into a single function call. Using the struct/functor approach, we can initialize each function with parameters, then pass the function and settings into the `floe_tracker` function. The default `FilterFunction` is an instance of a `ChainedFilterFunction` with seven inidividual filter functions, based on Lopez-Acosta et al. 2019.

After the filter function is applied, any floes remaining in ``C`` are given a column `trajectory_uuid` matching the `head_uuid` of the initial floe ``f``.

## Matching Function
Inevitably there will be cases where multiple candidates are plausible. `AbstractFloeMatchingFunctions` accept a DataFrame with the candidate pairings, then filter this dataframe to insure that all pairings are unique.

The default matching function is the `MinimumWeightMatchingFunction`. To initialize the function, it needs a list of properties. In this example, just three properties are assigned. The properties can be written as `Symbol` or as `String`.
```julia
matchfun = MinimumWeightMatchingFunction(
    columns=[:scaled_distance, :relative_error_area, :relative_error_convex_area]
)
```
Importantly, the names of the listed properties must coincide with columns created in the filter functions. They are interpreted as errors, such that a small value is favored over a large value.

The algorithm then works as follows. For each floe pair ``f, g``, and list of error measures ``m_1(f, g), \dots, m_k(f, g)``, we compute the weight
```math
w(f, g) = \sum_{i=1:k} m_i(f, g)
```
The matched pairs are then identified by finding the minimum ``w`` for each floe. We use a consistency check, such that we compute the minimum weight pair forward and backward and only return the matchings that minimize both directions.
