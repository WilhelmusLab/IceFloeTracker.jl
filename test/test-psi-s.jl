@testset "ψ-s curve test" begin

    # draw cardiod https://en.wikipedia.org/wiki/Cardioid
        t = range(0,2pi,201);
        x = @. cos(t)*(1-cos(t));
        y = @. (1-cos(t))*sin(t);

        ψ_s = IceFloeTracker.make_psi_s(x,y)
        
    # Test 0: Continuity for smooth curves (except initial/final point for cardiod)
        @test all((ψ_s[2:end] - ψ_s[1:end-1]) .< .05)

    # Test 1: initial/final phase [0, 3π]
        @test all([ψ_s[1]<0.05, abs(ψ_s[end]-3pi)<0.05])

    # Test 2: Option unwrap=false for the cardiod will have one discontinuity
        θ_s = IceFloeTracker.make_psi_s(x,y,rangeout=0,unwrap=false)
        @test sum((abs.(θ_s[2:end] - θ_s[1:end-1])) .> .05) == 1
    
    # Test 3: Catch bad rangeout values (good values for rangeout are 0,1)
        try
            IceFloeTracker.make_psi_s(x,y; rangeout=5)
        catch err
            @test isa(err, AssertionError)
        end

    # Test 4: Auxiliary functions
        A = [x y];
        # norm of a Vector
            @test all([IceFloeTracker.norm(x) == sum(x.^2)^.5, IceFloeTracker.norm(y) == sum(y.^2)^.5])

        # norm of a row/col
            @test all([IceFloeTracker.norm(A[1,:]) == sum(A[1,:].^2)^.5, IceFloeTracker.norm(A[:,1]) == sum(A[:,1].^2)^.5])    

        # test grad methods for vectors and matrices
            @test IceFloeTracker.grad(x,y) == IceFloeTracker.grad(A)

        # test arclength   
            @test sum([IceFloeTracker.norm(collect(v_i)) for v_i in eachrow(IceFloeTracker.grad(A))]) == IceFloeTracker.arclength(A)[end]     
end
