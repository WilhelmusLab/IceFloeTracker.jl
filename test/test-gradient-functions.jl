@testitem "gradient function – exponential" begin
    # Test basic functionality
    using IceFloeTracker: exponential
    @test exponential(0, 1) ≈ 1.0
    @test exponential(1, 1) ≈ exp(-1)
    @test exponential(2, 1) ≈ exp(-4)

    # Test with different k values
    @test exponential(1, 2) ≈ exp(-0.25)
    @test exponential(2, 2) ≈ exp(-1)

    # Test that function is monotonically decreasing
    @test exponential(0, 1) > exponential(1, 1)
    @test exponential(1, 1) > exponential(2, 1)

    # Test edge cases
    @test exponential(0, 1) == 1.0
    @test exponential(Inf, 1) ≈ 0.0 atol = 1e-10
end

@testitem "gradient function – inverse_quadratic" begin
    using IceFloeTracker: inverse_quadratic
    # Test basic functionality
    @test inverse_quadratic(0, 1) ≈ 1.0
    @test inverse_quadratic(1, 1) ≈ 0.5
    @test inverse_quadratic(2, 1) ≈ 0.2

    # Test with different k values
    @test inverse_quadratic(1, 2) ≈ 1 / (1 + 0.25)
    @test inverse_quadratic(2, 2) ≈ 1 / (1 + 1)

    # Test that function is monotonically decreasing
    @test inverse_quadratic(0, 1) > inverse_quadratic(1, 1)
    @test inverse_quadratic(1, 1) > inverse_quadratic(2, 1)

    # Test edge cases
    @test inverse_quadratic(0, 1) == 1.0
    @test inverse_quadratic(Inf, 1) ≈ 0.0 atol = 1e-10
end

@testitem "gradient function – properties" begin
    using IceFloeTracker: exponential, inverse_quadratic

    # Both functions should return values in (0, 1]
    test_values = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]
    k_values = [0.5, 1.0, 2.0, 5.0]

    for k in k_values, val in test_values
        exp_result = exponential(val, k)
        iq_result = inverse_quadratic(val, k)

        @test 0 < exp_result <= 1
        @test 0 < iq_result <= 1
    end
end

@testitem "gradient function – SupportedFunctions construction" begin
    using IceFloeTracker: SupportedFunctions

    sf = SupportedFunctions()

    @test sf isa SupportedFunctions
    @test length(sf.functions) == 2
    @test haskey(sf.functions, "exponential")
    @test haskey(sf.functions, "inverse_quadratic")
end

@testitem "gradient function – SupportedFunctions Base method overloads" begin
    using IceFloeTracker: SupportedFunctions

    sf = SupportedFunctions()

    # Test `in` operation
    @test "exponential" in sf
    @test "inverse_quadratic" in sf
    @test "nonexistent" ∉ sf

    # Test indexing
    @test sf["exponential"] isa Function
    @test sf["inverse_quadratic"] isa Function
    @test_throws KeyError sf["nonexistent"]

    # Test keys
    function_names = collect(keys(sf))
    @test "exponential" in function_names
    @test "inverse_quadratic" in function_names
    @test length(function_names) == 2
end

@testitem "gradient function – SupportedFunctions function retrieval and execution" begin
    using IceFloeTracker: SupportedFunctions, exponential, inverse_quadratic

    sf = SupportedFunctions()
    exp_func = sf["exponential"]
    iq_func = sf["inverse_quadratic"]

    # Test that retrieved functions work correctly
    @test exp_func(1, 1) ≈ exponential(1, 1)
    @test iq_func(1, 1) ≈ inverse_quadratic(1, 1)

    # Test with arrays
    test_array = [0.0, 1.0, 2.0]
    @test exp_func.(test_array, 1) ≈ exponential.(test_array, 1)
    @test iq_func.(test_array, 1) ≈ inverse_quadratic.(test_array, 1)
end

@testitem "gradient function – SupportedFunctions string representation" begin
    using IceFloeTracker: SupportedFunctions

    sf = SupportedFunctions()
    str_repr = string(sf)
    @test occursin("exponential", str_repr)
    @test occursin("inverse_quadratic", str_repr)
    @test occursin(",", str_repr)  # Should be comma-separated
end

@testitem "gradient function – is_supported helper function" begin
    using IceFloeTracker: SupportedFunctions, is_supported

    sf = SupportedFunctions()
    @test is_supported(sf, "exponential")
    @test is_supported(sf, "inverse_quadratic")
    @test !is_supported(sf, "nonexistent")
    @test !is_supported(sf, "")
end

@testitem "gradient function – SUPPORTED_GRADIENT_FUNCTIONS constant" begin
    using IceFloeTracker: SUPPORTED_GRADIENT_FUNCTIONS, SupportedFunctions

    @test SUPPORTED_GRADIENT_FUNCTIONS isa SupportedFunctions
    @test "exponential" in SUPPORTED_GRADIENT_FUNCTIONS
    @test "inverse_quadratic" in SUPPORTED_GRADIENT_FUNCTIONS

    # Test that it behaves the same as a new instance
    new_sf = SupportedFunctions()
    @test collect(keys(SUPPORTED_GRADIENT_FUNCTIONS)) == collect(keys(new_sf))
end

@testitem "gradient function – Integration with nonlinear diffusion" begin
    using IceFloeTracker: SUPPORTED_GRADIENT_FUNCTIONS

    # Test that the functions work as expected in the context they'll be used
    # Simulate typical values that might come from image gradients
    gradient_magnitudes = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    k_values = [1, 5, 10, 50, 100]

    for k in k_values
        exp_func = SUPPORTED_GRADIENT_FUNCTIONS["exponential"]
        iq_func = SUPPORTED_GRADIENT_FUNCTIONS["inverse_quadratic"]

        exp_results = exp_func.(gradient_magnitudes, k)
        iq_results = iq_func.(gradient_magnitudes, k)

        # Results should be decreasing as gradient magnitude increases
        @test issorted(exp_results, rev=true)
        @test issorted(iq_results, rev=true)

        # All results should be in (0, 1]
        @test all(0 < x <= 1 for x in exp_results)
        @test all(0 < x <= 1 for x in iq_results)
    end
end
