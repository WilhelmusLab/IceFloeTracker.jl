@testitem "ψ-s curve test" begin
    # draw cardioid https://en.wikipedia.org/wiki/Cardioid
    t = range(0, 2pi, 201)
    x = @. cos(t) * (1 - cos(t))
    y = @. (1 - cos(t)) * sin(t)

    ψ, s = buildψs(x, y)

    # Test 0: Continuity for smooth curves (except initial/final point for cardioid)
    @test all((ψ[2:end] - ψ[1:(end - 1)]) .< 0.05)

    # Test 1: initial/final phase [0, 3π]
    @test all([ψ[1] < 0.05, abs(ψ[end] - 3pi) < 0.05])

    # Test 2: Option unwrap=false; cardioid will have one discontinuity
    θ, s = buildψs(x, y; rangeout=true, dsp_unwrap=false)
    @test sum((abs.(θ[2:end] - θ[1:(end - 1)])) .> 0.05) == 1

    # Test 3: Option rangeout=false, unwrap=false. There will be negative phase values.
    p_wrapped, _ = buildψs(x, y; rangeout=false, dsp_unwrap=false)
    @test !all(p_wrapped .>= 0)

    # Test 4: Return arclength; compare against theoretical total arclength = 8
    _, s = buildψs(x, y; rangeout=true)
    @test abs(s[end] - 8) < 0.005

    # Test 5: Auxiliary functions
    A = [x y]
    # norm of a Vector
    @test all([isapprox(norm(x), sum(x .^ 2)^0.5), isapprox(norm(y), sum(y .^ 2)^0.5)])

    # norm of a row/col
    @test all([
        isapprox(norm(A[1, :]), sum(A[1, :] .^ 2)^0.5),
        isapprox(norm(A[:, 1]), sum(A[:, 1] .^ 2)^0.5),
    ])

    # test grad methods for vectors and matrices
    @test grad(x, y) == grad(A)

    # Test 6: Alternate method of make_psi_s with 2-column matrix as input
    @test buildψs(x, y) == buildψs([x y])
end
