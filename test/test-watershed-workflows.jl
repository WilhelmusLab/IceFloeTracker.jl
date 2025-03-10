@testset "watershed workflows" begin
    function build_test_image()
        center1, center2 = -40, 40
        radius = sqrt(8 * center1^2) * 0.7
        lims = (floor(center1 - 1.2 * radius), ceil(center2 + 1.2 * radius))
        x = collect(lims[1]:lims[2])

        function in_circle(x, y, center, radius)
            return sqrt((x - center)^2 + (y - center)^2) <= radius
        end

        return [
            in_circle(xi, yi, center1, radius) || in_circle(xi, yi, center2, radius) for
            yi in x, xi in x
        ]
    end
    @testset "watershed" begin
        println("------------------------------------------------")
        println("------------ Create Watershed Test --------------")
        @test sum(IceFloeTracker.watershed1(build_test_image())) == 1088
    end

    # TODO: #589 add test for watershed2
end
