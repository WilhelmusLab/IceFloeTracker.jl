@testset "bwperim" begin
    println("-------------------------------------------------")
    println("---------------- bwperim Tests ------------------")

    # Create image with 3 connected components. The test consists of digging the biggests holes for each blob in the foreground using bwperim, thereby creating 3 additional connected components, 6 in total.
    A = zeros(Bool, 13, 16)
    A[2:6, 2:6] .= 1
    A[4:8, 7:10] .= 1
    A[10:12, 13:15] .= 1
    A[10:12, 3:6] .= 1
    # 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    # 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0
    # 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0
    # 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
    # 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
    # 0  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0
    # 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0
    # 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0
    # 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    # 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
    # 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
    # 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
    # 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

    # count the components in A
    numinicomps = maximum(IceFloeTracker.label_components(A))

    # dig the holes and relabel the interiors
    border_mask = IceFloeTracker.bwperim(A) # get the borders
    interior = A .& .!border_mask # keep the interior of the holes
    interior_relabel = IceFloeTracker.label_components(interior) * 2 # relabel the wholes
    A_relabeled = interior_relabel .+ border_mask
    # 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    # 0  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0        
    # 0  1  2  2  2  1  0  0  0  0  0  0  0  0  0  0        
    # 0  1  2  2  2  2  1  1  1  1  0  0  0  0  0  0        
    # 0  1  2  2  2  2  2  2  2  1  0  0  0  0  0  0        
    # 0  1  1  1  1  1  2  2  2  1  0  0  0  0  0  0        
    # 0  0  0  0  0  0  1  2  2  1  0  0  0  0  0  0        
    # 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0        
    # 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0        
    # 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
    # 0  0  1  4  4  1  0  0  0  0  0  0  1  6  1  0        
    # 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0        
    # 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

    @test numinicomps + 3 ==
          maximum(IceFloeTracker.label_components(A_relabeled, trues(3, 3)))
end
