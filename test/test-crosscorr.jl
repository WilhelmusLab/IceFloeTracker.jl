@testset "crosscorr tests" begin
    println("-------------------------------------------------")
    println("----------- cross correlation tests -------------")

    # Compare against matlab xcorr with normalization (standard use case in IceFloeTracker)
        n = 0:15;
        x = 0.84.^n;
        y = circshift(x,5);
        # [cnormalized,lags] = xcorr(x,y,'coef'); (matlab call)
        
    # Output of matlab call above
        c_matlab = [0.0516860456455592,0.104947285063174,0.161406917930332,0.222785618772511,0.290953976567756,0.367989503172686,0.456239947969545,0.558394848323571,0.677567496436029,0.817389820630348,0.982123072691496,0.846599102605261,0.736876248026996,0.649610574340982,0.582142556253932,0.532416025595572,0.444053743842906,0.369224528569261,0.305647870356775,0.251386194859924,0.204785812920709,0.164426522422887,0.129078325941762,0.097663945108386,0.0692259892687896,0.0428977778640515,0.0178769273085035,0.0136884958891382,0.00991723767782281,0.00644821909097442,0.00317571765737478]

    # Compute normalized cross correlation scores and lags with padding
        r, lags = IceFloeTracker.crosscorr(x,y,normalize=true)

    # Test
        @test all(round.(c_matlab,digits=5) .== round.(r,digits=5))
end