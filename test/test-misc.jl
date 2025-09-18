# Define a test to check the version number of the package matches the version number in the Project.toml file. Use the get_version_from_toml function to get the version number from the Project.toml file. The version number is then compared to the version number of the package. If the version numbers match, the test passes. If the version numbers do not match, the test fails.

@testitem "Version number" begin
    VERSION = IceFloeTracker.get_version_from_toml(dirname((@__DIR__)))
    @test isa(VERSION, VersionNumber)
end
