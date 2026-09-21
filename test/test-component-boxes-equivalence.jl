@testitem "_component_boxes agrees with ImageMorphology" begin
    using IceFloeTracker.Segmentation: _component_boxes
    using Images: component_boxes

    # TODO: Remove this test once the fixed version of `ImageMorphology.component_boxes`
    # is released and the local copy is no longer needed.
    #
    # src/Segmentation/component-boxes.jl is a local copy of `ImageMorphology.component_boxes` carrying the
    # fix in JuliaImages/ImageMorphology.jl#146, which changes that function's
    # speed and not its results, which is asserted here. Once a release containing
    # the fix is available, both go: the copy and these tests.


    A = [2 2 2 2 2; 1 1 1 0 1; 1 0 2 1 1; 1 1 2 2 2; 1 0 2 2 2] # from the function's docstring
    ours, theirs = _component_boxes(A), component_boxes(A)
    @test axes(ours) == axes(theirs) && all(ours[i] == theirs[i] for i in eachindex(ours))

end
