using Images
# @testset "pipeline" begin
    pipelinedir = joinpath(test_data_dir,"pipeline")
    input = joinpath(pipelinedir,"input")
    output = joinpath(pipelinedir, "output")
    args = Dict(zip([:input, :output],[input, output]))
    IceFloeTracker.landmask(args...)
    lm_expected = Gray.(load(joinpath(test_data_dir,"pipeline","expected","generated_landmask.png"))).>0

    # end